//https:www.cairographics.org/manual/cairo-Paths.html fix
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>
#include "Graphics.h"

void
Graphics::draw(Particle *p, double FPS,int l){
  
  //  const int FPS=50;//frames per second
  struct timespec tm={0,1000000000/FPS}; //Sleep in nanoseconds between frames
  XEvent event;//check if window closed and finish
  if( XCheckWindowEvent(dsp,da, DestroyNotify , &event)){ XCloseDisplay(dsp); exit(1);}

  double gamma = Dim/(lmax-lmin + diam);//scaling between physical units and pixels
  double alpha= gamma*diam/2.;
  char a[10];
      
  cairo_push_group(cr); //start drawing
  cairo_set_source_rgb(cr, 0.89, 0.89, .89);//dark green background
  cairo_paint (cr); //clear screen with green

  if(1){//solid particles
    //dessiner la paroi 
    cairo_set_source_rgb(cr, 0.25, 0.25, 0.25);
    cairo_set_line_width(cr, 5);
    cairo_move_to(cr,(alpha + gamma* (p[Np].x -lmin)),Dim-(alpha + gamma*(0 - lmin)));
    cairo_line_to(cr,(alpha + gamma* (p[Np].x -lmin)),Dim-(alpha + gamma*(lmax-lmin)));
    cairo_stroke(cr);

    cairo_set_source_rgb(cr, 0.58, 0.065, 0.52);//dark red for particles  
    for(int i=0;i<Np;i++){// place the particles in the graphics buffer, without drawing
      cairo_new_sub_path(cr);
      cairo_arc(cr,   (alpha + gamma* (p[i].x -lmin)) ,  Dim-(alpha + gamma*(p[i].y - lmin)), alpha, 0, 2 * M_PI);
    }
      cairo_fill(cr);//draw all particles with solid color
  }
  
  
  if(0){//label particles with numbers
    cairo_set_source_rgb(cr, 1, 1, 1);//white
    for(int i=0;i<Np;i++){// number particles
  
      cairo_move_to (cr, (alpha + gamma* (p[i].x -lmin)) ,  Dim-(alpha + gamma*(p[i].y - lmin)));
      sprintf(a, "%d",i);
      cairo_set_font_size(cr,15);

      cairo_show_text (cr, a);
    }
  }
  
  cairo_move_to(cr, Dim-120, Dim-35);//print frame number to graphics window
  sprintf(a, "%d",l);
  cairo_set_font_size(cr,20);
  cairo_show_text (cr, a);

  if(0){//RING for outter circle
    cairo_set_source_rgb(cr, .2, 0.79, 0.79);//dark blue for particles
    for(int i=0;i<Np;i++){
      cairo_new_sub_path(cr) ;
      cairo_arc(cr,  (alpha + gamma* (p[i].x -lmin)) ,  Dim - (alpha + gamma*(p[i].y - lmin)), 1.4*alpha, 0, 2 * M_PI);
    }
    cairo_stroke (cr); // hollow particles
  }
  
  cairo_set_line_width(cr,1);
  frame(alpha, Dim-alpha, alpha, Dim-alpha);//draw square border
  cairo_pop_group_to_source(cr); //finished drawing operations for this set of positions
  cairo_paint(cr);//send to screen
  cairo_surface_flush(sfc); //send to x11
  nanosleep( &tm , NULL); //this sets the animations speed
  XFlush(dsp);//sync X11 to cairo


  //https://cairographics.org/Xlib/

  
}

void Graphics::frame(double xmin, double xmax, double ymin, double ymax){
  cairo_set_source_rgb (cr, 0, 0, 0);//yellow border
  cairo_rectangle (cr,xmin, ymin, 3*(xmax-xmin), ymax-ymin);
  cairo_stroke (cr);
}
