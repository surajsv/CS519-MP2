
//University of Illinois/NCSA Open Source License
//Copyright (c) 2015 University of Illinois
//All rights reserved.
//
//Developed by: 		Eric Shaffer
//                  Department of Computer Science
//                  University of Illinois at Urbana Champaign
// Edited and modified by: Suraj Venkat
//                          
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
//documentation files (the "Software"), to deal with the Software without restriction, including without limitation
//the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and
//to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
//Redistributions of source code must retain the above copyright notice, this list of conditions and the following
//disclaimers.Redistributions in binary form must reproduce the above copyright notice, this list
//of conditions and the following disclaimers in the documentation and/or other materials provided with the distribution.
//Neither the names of <Name of Development Group, Name of Institution>, nor the names of its contributors may be
//used to endorse or promote products derived from this Software without specific prior written permission.
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
//WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
//TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//DEALINGS WITH THE SOFTWARE.




//-------------------------------------------------------
// Global variables

var x_extent=[-1.0,1.0];
var y_extent=[-1.0,1.0];
var noise_arr = [] ;
var LIC = [] ;
var mag_arr = [] ;
var myGrid;
var flag = 0;
var color_switch = 0;
var line_length = 10;
var k_factor = 1;
var grid_sz = 10;
var hedgehog_switch = 0;


function gaussian(pt){
	return Math.exp(-(pt[0]*pt[0]+pt[1]*pt[1]));
}
function weight(pt1,pt2)
{
    var disyy = Math.abs(pt1[1]-pt2[1]) ;
    var disxx = Math.abs(pt1[0]-pt2[0]) ;
    return Math.exp(-(disyy*disyy+disxx*disxx)) ;
    
}
//--------------------------------------------------------
//The infamous rainbow color map, normalized to the data range
function rainbow_colormap(fval,fmin,fmax){
	var dx=0.8;
	var fval_nrm = (fval-fmin)/(fmax-fmin);
	var g = (6.0-2.0*dx)*fval_nrm +dx;
	var R = Math.max(0.0,(3.0-Math.abs(g-4.0)-Math.abs(g-5.0))/2.0 )*255;
	var G = Math.max(0.0,(4.0-Math.abs(g-2.0)-Math.abs(g-4.0))/2.0 )*255;
	var B = Math.max(0.0,(3.0-Math.abs(g-1.0)-Math.abs(g-2.0))/2.0 )*255;
	color = [Math.round(R),Math.round(G),Math.round(B),255];
	return color;
}

//--------------------------------------------------------
//

function greyscale_map(fval,fmin,fmax){
  var c=255*((fval-fmin)/(fmax-fmin));
  var color = [Math.round(c),Math.round(c),Math.round(c),255];
	return color;
}

//--------------------------------------------------------
//

function gaussian_gradient(pt){
  var dx = -2*pt[0]*gaussian(pt);
  var dy = -2*pt[1]*gaussian(pt);
    
    var ret = normalize2D([dx,dy]) ;
	return (ret);
}

function gaussian_gradient_mag(pt){
  var dx = -2*pt[0]*gaussian(pt);
  var dy = -2*pt[1]*gaussian(pt);
    var mag = (dx*dx) + (dy*dy) ;
    var ret = Math.sqrt(mag) ;
	return (ret);
}

//--------------------------------------------------------
//

function gaussian_divergence(pt){
  var gradient =  gaussian_gradient(pt);
	return gradient[0]+gradient[1];
}

//--------------------------------------------------------
//

function gaussian_vorticity_mag(pt){
	return 0;
}

//--------------------------------------------------------
//

function normalize2D(v){
  var len = Math.sqrt(v[0]*v[0] + v[1]*v[1]);
  if (len == 0.0)
    {
       console.log("Zero length gradient");
       return ([0.0,0.0]);
    }
  return [v[0]/len,v[1]/len];
}

//--------------------------------------------------------
//

function euler_integration(pt,h,steps,get_vector)
{
    var ln=[[pt[0],pt[1]]];
    for(i=0;i<steps;i++)
      {
        v = get_vector(ln[i]);
        ln.push( [ln[i][0]+h*v[0],ln[i][1]+h*v[1]]);
      }
    return ln;
}

//--------------------------------------------------------
// Map a point in pixel coordinates to the 2D function domain
function pixel2pt(width,height,x_extent,y_extent, p_x,p_y){
	var pt = [0,0];
	xlen=x_extent[1]-x_extent[0]
	ylen=y_extent[1]-y_extent[0]
	pt[0]=(p_x/width)*xlen + x_extent[0];
	pt[1]=(p_y/height)*ylen + y_extent[0];
	return pt;
	}

//--------------------------------------------------------
// Map a point in domain coordinates to pixel coordinates
function pt2pixel(width,height,x_extent,y_extent, p_x,p_y){
	var pt = [0,0];

	var xlen = (p_x-x_extent[0])/(x_extent[1]-x_extent[0]);
  var ylen = (p_y-y_extent[0])/(y_extent[1]-y_extent[0]);

	pt[0]=Math.round(xlen*width);
	pt[1]=Math.round(ylen*height);
	return pt;
	}


//------------------------------------------------------
//MAIN
function main() {
	render();
}

//--Function: render-------------------------------------
//Main drawing function
var UGrid2D = function(min_corner,max_corner,resolution){
  this.min_corner=min_corner;
  this.max_corner=max_corner;
  this.resolution=resolution;
  console.log('UGrid2D instance created');
}


// Method: draw_grid
// Draw the grid lines

UGrid2D.prototype.draw_grid = function(canvas){
	  var ctx = canvas.getContext('2d');
      
//This is to draw the LIC

var L = line_length; //Length
  
var lim = Math.ceil(L/2) ;      
for (var y=0;y<canvas.height;y++)
{

    for (var x=0;x<canvas.width;x++)
    {
        var source_point = pixel2pt(canvas.width,canvas.height,x_extent,y_extent,x,y); //this is the source point
        
        var sum_d=0;
        var sum_n = 0;  // first two for summing numerator forward streamline
  		
        var sum_d2=0;  //next two for backward streamline
        var sum_n2 = 0;
        
        sum_d = 1;
        sum_n = noise_arr[(y*canvas.width) + x] ;
        
        sumd_d2 = 1;
        sum_n2= noise_arr[(y*canvas.width) + x];
        var weight_func= 0;
        
        var dx = gaussian_gradient(source_point)[0] ;
        var dy = gaussian_gradient(source_point)[1] ;
        
       
        
        var new_x = source_point[0] ;  //these move in forward stream
        var new_y = source_point[1] ;
        
        var new_x2 = source_point[0] ; //these move in backward stream
        var new_y2 = source_point[1] ;
        
        var h = 2/3; //step size
        var x_lim = canvas.width;
        x_lim = x_lim - lim;
        var y_lim = canvas.height;
        y_lim = y_lim - lim;
     
        
        
        for(var i=1; i<=L;i++)
        {
        
            dx = gaussian_gradient([new_x,new_y])[0] * (2/canvas.width) ;//dx and dy for forward stream (gaussian_gradient handles the conversio 
                                                        //-n to unit vector
            
            dy = gaussian_gradient([new_x,new_y])[1] * (2/canvas.width) ;
            
            neg_dx = - (gaussian_gradient([new_x2,new_y2])[0]) * (2/canvas.width)  ; //these are for moving in backward direction
            neg_dy = -(gaussian_gradient([new_x2,new_y2])[1]) * (2/canvas.width) ;
                    
            
            new_x = new_x + (h*dx); //to move by step size
          
            new_y = new_y + (h*dy);
                
            new_x2 = new_x2 + (h*neg_dx);
                
            new_y2 = new_y2 + (h*neg_dy) ;    
                
            
      
            //console.info([new_x,new_y])
            if(new_x > 1.0 || new_y > 1.0 || new_x< -1.0 || new_y < -1.0 )
            {    
              sum_n = sum_n ;
              sum_d = sum_d ; 
            }
            
             else{
                var new_pix = pt2pixel(canvas.width,canvas.height,x_extent,y_extent,new_x,new_y);
                
                if((0<new_pix[0]<canvas.width) && (0<new_pix[1]<canvas.height))
                {
                
                var disx = new_x-source_point[0];
                var disy = new_y-source_point[1] ;
                
                weight_func = gaussian([disx,disy]);
                
                var ind = (new_pix[1]*canvas.width) + new_pix[0];
                if(ind<160000)
                {
                if(ind >=0)
                {
                sum_n = sum_n + (weight_func * noise_arr[ind]);
          
                sum_d  = sum_d + weight_func ;
                }
                }
                }
            }
            
            if(new_x2 > 1.0 || new_y2 > 1.0 || new_x2 < -1.0 || new_y2 < -1.0 )
            {    
              sum_n2 = sum_n2 ;
              sum_d2 = sum_d2 ; 
            }
            
             else{
                var new_pix = pt2pixel(canvas.width,canvas.height,x_extent,y_extent,new_x2,new_y2);
                
                 
                if((0<new_pix[0]<canvas.width) && (0<new_pix[1]<canvas.height))
                {
                
                var disx = new_x2-source_point[0];  //distance in terms of points and not pixels
                var disy = new_y2-source_point[1] ;
                weight_func = gaussian([disx,disy]) ;
                
                
                var ind = (new_pix[1]*canvas.width) + new_pix[0];
                if(ind<160000)
                {
                    if(ind>=0)
                    {
                    sum_n2 = sum_n2 + (weight_func * noise_arr[ind]);
          
                    sum_d2  = sum_d2 + weight_func ;
                    }
                }
                }
          
            }
           
        }
       
        var num = sum_n + sum_n2 ; //summing numerator terms
        var den = sum_d + sum_d2 ; //summing denominator terms
        var ans = 0;
        if(den!=0)
         ans =  num/den ; //this is the LIC for that pixel!
        else
         ans = 0;
        LIC[y*canvas.width + x] = ans ; //row major indexing
    }
    
console.log("one set done");
}
    console.log("noise")
    console.info(noise_arr);
    console.log("LIC:")
    console.info(LIC);
    
var imgData=ctx.getImageData(0,0,canvas.width,canvas.height); 
    

    
var min = LIC[0];
var max = LIC[0];
var pt1 = [-1.0,-1.0];
    
for (var y=0;y<canvas.height;y++)
{
    for (var x=0;x<canvas.width;x++)
  	{
  		var f = LIC[y*canvas.width + x];
        pt1 = pixel2pt(canvas.width,canvas.height,x_extent,y_extent,x,y);
  		var mag = gaussian_gradient_mag(pt1);
        
        mag_arr.push(mag) ;
        
        if (f < min)
  			min=f;
  		if (f>max)
  			max=f;
  	}
}
    
    
    
var min2 = mag_arr[0];
var max2 = mag_arr[0];


    
for (var y=0;y<canvas.height;y++)
{
	for (var x=0;x<canvas.width;x++)
  	{
  		var f = mag_arr[y*canvas.width + x];
        if (f < min2)
  			min2=f;
  		if (f>max2)
  			max2=f;    
    }    

}    
    console.log("min2:");
    console.info(min2);
    console.log("max2:");
    console.info(max2);
    
for (var y=0;y<canvas.height;y++)
{
	for (var x=0;x<canvas.width;x++)
  	{
  		var fff = LIC[ y*canvas.width + x ]; //extract the LIC value
        var f = mag_arr[ y*canvas.width + x ] ;
        if(color_switch==1)
        {
        var color = rainbow_colormap(f,min2,max2); //now map to gray scale
        i = (y*canvas.width +x ) * 4 ;
  		var temp = color[0] * fff ;
        imgData.data[i]= temp ;
        temp = color[1] * fff ;
  		imgData.data[i+1]= temp ;
        temp = color[2] * fff;
  		imgData.data[i+2]= temp ;
        temp = color[3] * fff ;
  		imgData.data[i+3]= temp ;    
  		
        }
        else
        {
        var color2 = grey_func(fff,min,max);
        i = (y*canvas.width +x ) * 4 ;
  		temp = color2[0];
        imgData.data[i]= temp ;
        temp = color2[1] ;
  		imgData.data[i+1]= temp ;
        temp = color2[2];
  		imgData.data[i+2]= temp ;
        temp = color2[3] ;
  		imgData.data[i+3]= temp ;    
        }
    }    

}    
       
    

    ctx.putImageData(imgData,0,0);    
    
    
    //The following portion is for hedgehog plot which works
    if(hedgehog_switch==1)
    {
     loc=[0,0];
    var res = parseFloat(document.getElementById("grid_size").value);
	  var deltax = canvas.width/(res);
      var deltay =canvas.height/(res);

    for( var j=0; j< res; j++)
	  for (var i=0;i< res; i++)
	  {
      ctx.beginPath();
        var startx = Math.floor((2*i+1)*deltax*0.5);
        var starty = Math.floor((2*j+1)*deltay*0.5);
          console.log("start:");
        console.info([startx,starty]);
	  	ctx.moveTo(startx, starty);
        
        var pt = pixel2pt(canvas.width,canvas.height,x_extent,y_extent,startx,starty);
          
        var end = gaussian_gradient(pt) ;
        end[0] = pt[0] + k_factor*10*(end[0] * (2/canvas.width)) ;
        end[1] = pt[1] + k_factor*10*(end[1] * (2/canvas.width)) ;
          
        end = pt2pixel(canvas.width,canvas.height,x_extent,y_extent,end[0],end[1]);
        console.log("end:");
        console.info(end);
        ctx.lineTo(end[0], end[1]);
      	ctx.lineWidth = 2;
      	// set line color
      	ctx.strokeStyle = '#FF0000';
      	ctx.stroke();
	   }   
    }
}




//End UGrid2D--------------------------------------------
function render(canvas){
  var res = parseFloat(document.getElementById("grid_size").value);
  line_length = parseFloat(document.getElementById("line_length").value);
  k_factor = parseFloat(document.getElementById("k_hedgehog").value);
  myGrid = new UGrid2D([x_extent[0],y_extent[0]],  [x_extent[1],y_extent[1]]  ,res);
  var canvas = document.getElementById('example');
  if (! canvas) {
    console.log(' Failed to retrieve the < canvas > element');
    return false;
  }
  else {
	console.log(' Got < canvas > element ');
  }
var ctx = canvas.getContext('2d');
if(flag==1)
myGrid.draw_grid(canvas);
// Get the rendering context for 2DCG <- (2)

//var cty = canvas.getContext('2d');
// Draw the scalar data using an image rpresentation
var imgData=ctx.getImageData(0,0,canvas.width,canvas.height);
//var imgData2=cty.getImageData(0,0,canvas.width,canvas.height);

if (document.getElementById("hedge_on").checked)
            hedgehog_switch = 1;

else
    hedgehog_switch = 0;
    
if (document.getElementById("color_on").checked)
            color_switch = 1;

else
    color_switch = 0; 
    
// Choose the scalar function
var scalar_func = gaussian;
var reso = 2/(canvas.width);
reso = reso -1 ;

//Determine the data range...useful for the color mapping
var mn = scalar_func(pixel2pt(canvas.width,canvas.height,x_extent,y_extent,0,0));
var mx = mn ;
var noise_min = 0;
var noise_max = 1.0;
if(flag==0)
{
flag = 1;
for (var y=0;y<canvas.height;y++)
{
	for (var x=0;x<canvas.width;x++)
  	{
  		var fval = scalar_func(pixel2pt(canvas.width,canvas.height,x_extent,y_extent,x,y));
  		if (fval < mn)
  			mn=fval;
  		if (fval>mx)
  			mx=fval;
        var noise = Math.random();
        noise_arr.push(noise);
        LIC.push(noise);
  	}
}
console.log("noise:");
    console.info(noise_arr);
mn = noise_arr[0];
mx = noise_arr[0];
    
for (var y=0;y<canvas.height;y++)
{
    for (var x=0;x<canvas.width;x++)
  	{
  		var fval = noise_arr[y*canvas.width + x];
  		if (fval < mn)
  			mn=fval;
  		if (fval>mx)
  			mx=fval;
        
  	}
}
console.log("min:");
console.info(mn);
    
console.log("max:");
console.info(mx);
    
// Set the colormap based in the radio button
var color_func = rainbow_colormap;


grey_func = greyscale_map;
    
//Color the domain according to the scalar value
for (var y=0;y<canvas.height;y++)
{
    for (var x=0;x<canvas.width;x++)
  	{
  		var fval = scalar_func(pixel2pt(canvas.width,canvas.height,x_extent,y_extent,x,y));
        
        
  		//var color = color_func(fval,mn,mx);
        var noisy = grey_func(noise_arr[y*canvas.width +x],mn,mx);
    
        i = (y*canvas.width +x ) * 4 ;
  		
        imgData.data[i]=noisy[0];
  		imgData.data[i+1]= noisy[1];
  		imgData.data[i+2]= noisy[2];
  		imgData.data[i+3]= noisy[3]; 
        
    }
}


ctx.putImageData(imgData,0,0);
    
}
//NOW DO LIC calculation


  // Draw the grid if necessary
  
    
}


