#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>

#include "Graphics.h"  // Contains the definition of "Particle" and "Graphics"

enum  col_type {
  bottom,
  right,
  top,
  left,
  animation,
  particle,
  paroi,

}; //different types of collision


typedef struct { // a structure describing each collision -- one might want an array of Events
  enum col_type type;
  int ia;
  int ib;
  double time;
} Event;

typedef struct{ //structure permettant de memoriser le type de choc, notamment le choc avec la paroi
  enum col_type type;
  int ia;
  int ib;
} Mem_paroi;


void compute_parroi_collision_time(Particle *p, int Np, Event *e, double diameter){
  int compteur = 4*Np+1+(Np*(Np-1))/2;
  int i;

  for(i=0; i<Np;i++){ //va jusqu'a Np car pas de collision de la paroi avec elle meme
    e[compteur].ia = i;
    e[compteur].ib = Np;
    e[compteur].type = paroi;
    if(p[Np].x > p[i].x && (p[i].vx-p[Np].vx)>0){
      e[compteur].time = (p[Np].x - diameter/2 - p[i].x)/(p[i].vx-p[Np].vx);
    }

    else if( p[Np].x < p[i].x && (p[Np].vx-p[i].vx)>0){
      e[compteur].time = (p[Np].x + diameter/2 - p[i].x)/(p[i].vx-p[Np].vx);
    }

    compteur++;
  }//compteur est a 4*Np+1+(Np*(Np-1))/2 + Np
}


void compute_particle_collision_time(Particle *p, int Np, Event *e, double diameter){
  int compteur = 4*Np+1;
  int i, j;
  double a , b, c, DELTA;
  double dx , dy, dvx, dvy ;

  for (i = 0; i < Np; i++)
  {
    for (j = i+1; j < Np; j++){
      dx = p[i].x-p[j].x;
      dy = p[i].y-p[j].y;
      dvx = p[i].vx-p[j].vx;
      dvy = p[i].vy-p[j].vy;

      b = 2*( dx * dvx + dy * dvy);
      if( b < 0){ //verifie dans un premier temps si b<0 avant de faire le calcul.
        a = dvx * dvx + dvy * dvy;
        c = dx * dx + dy * dy - diameter * diameter;
        DELTA = b * b - 4. * a * c;
        if (DELTA > 0){
          e[compteur].time = (-b-sqrt(DELTA))/(2*a);
          e[compteur].type = particle;
          e[compteur].ia = i;
          e[compteur].ib = j;
        }
        else
        {
          e[compteur].time = -1;
        }
        
      }
      else{
        e[compteur].time = -1;
      }

      compteur++;
    }
  }
  //if faut compteur == 4*Np+(N*(N-1)/2
} 


void compute_wall_collision_time(Particle *p, int Np, Event *e, double diameter, double Lmax, double Lmin, double lgmax){ 
  int compteur=0;
  int i;
  for (i = 0; i < Np; i++){
    e[compteur].ia = i;
    e[compteur].type = top;
    if(p[i].vy>0){
      e[compteur].time = (Lmax - diameter/2 - p[i].y)/p[i].vy;
    }
    else{
      e[compteur].time = -1;
    }   
    compteur++; 
    e[compteur].ia = i;
    e[compteur].type = right;
    if(p[i].vx>0){
      e[compteur].time = (lgmax - diameter/2 - p[i].x)/p[i].vx;
    }
    else{
      e[compteur].time = -1;
    }   
    compteur++;
    e[compteur].ia = i;
    e[compteur].type = left;

    if(p[i].vx<0){
      e[compteur].time = (Lmin + diameter/2 - p[i].x)/p[i].vx;
    }
    else{
      e[compteur].time = -1;
    }   
    compteur++;  
    e[compteur].ia = i;
    e[compteur].type = bottom;
    if(p[i].vy<0){
      e[compteur].time = (Lmin + diameter/2 - p[i].y)/p[i].vy;
    }
    else{
      e[compteur].time = -1;
    } 
    compteur++;  
  }
}

int compute_min_time(int Np, Event *e, Mem_paroi M){
  int i, pos;
  double tmin = 10000000;
  
  //si une particule a effectué une collision avec la paroi au temps t, elle ne peut pas retoucher la paroi au temps t+1
  for (i = 0; i < (4*Np+1+(Np*(Np-1))/2+Np); i++){
    if(e[i].type == M.type && e[i].ia == M.ia && e[i].ib== M.ib && e[i].type==paroi){
      e[i].time = -1;
    } 
  }
  
  //calcul du min time
  for (i = 0; i < (4*Np+1+(Np*(Np-1))/2+Np); i++){
    if (e[i].time < tmin && e[i].time > 0){
      tmin = e[i].time;
      pos = i;
    }
  }
  return pos;  
}

void new_pos(Particle *p, int Np, double tmin){ //changer void en double pour pouvoir recuperer la position moyenne de la paroi
  for (int k = 0; k < Np; k++){ //on inclut la parroi
    p[k].x += p[k].vx * tmin;
    p[k].y += p[k].vy * tmin;
  }
  p[Np].x += p[Np].vx * tmin;
  p[Np].y = 0; //fixe parroi

//return p[Np].x; //decommenter pour recuperer position moyenne de la paroi
}

void new_speed(Particle *p, int pos, Event *e, int Np){
  
  if(e[pos].type == paroi){ //mise a jour des vitesses après collision avec paroi
    int particule = e[pos].ia;
    int mur = Np;
    double temp;
    temp = p[mur].vx;
    p[mur].vx = p[particule].vx ;
    p[particule].vx = temp;
    
  }
  
  if(e[pos].type == particle){ //mise a jour des vitesses après collision entre particules
    int particule_1 = e[pos].ia;
    int particule_2 = e[pos].ib;
    int nb_chocs = 0;

    double dx = (p[particule_2].x-p[particule_1].x);
    double dy = (p[particule_2].y-p[particule_1].y);
    double normed=sqrt(dx*dx+dy*dy);
    dx=dx/normed;
    dy=dy/normed;
    
    //projections
    double proj1x = (p[particule_1].vx*dx + p[particule_1].vy*dy)*dx;
    double proj1y = (p[particule_1].vx*dx + p[particule_1].vy*dy)*dy;

    double proj2x = (p[particule_2].vx*dx + p[particule_2].vy*dy)*dx;
    double proj2y = (p[particule_2].vx*dx + p[particule_2].vy*dy)*dy;

    double projperp1x = p[particule_1].vx - proj1x;
    double projperp1y = p[particule_1].vy - proj1y;

    double projperp2x = p[particule_2].vx - proj2x;
    double projperp2y = p[particule_2].vy - proj2y;

    p[particule_1].vx=projperp1x+proj2x;
    p[particule_1].vy=projperp1y+proj2y;

    p[particule_2].vx=projperp2x+proj1x;
    p[particule_2].vy=projperp2y+proj1y;

  }
    /*---------------------*/
    double norme_max = 5;     //faire varier pour changer le gradient, norme maximale de la vitesse qu'une particule peut prendre
    double norme_min = 0.3;  
    /*---------------------*/

    if(e[pos].type == left){ //mise a jour des vitesses avec mur "chaud"
      int particule = e[pos].ia;

      double norme = sqrt(pow(p[particule].vx,2)+pow(p[particule].vy,2));
      double rd = drand48(); //permet d'avoir le même drand48() pour vx et vy et ainsi conserver l'angle

      //on norme la vitesse
      p[particule].vx = p[particule].vx/norme; 
      p[particule].vy = p[particule].vy/norme;

      //on modelise maxwell-boltzmann par une loi de proba uniforme
      p[particule].vx *= -rd*norme_max;
      p[particule].vy *= rd*norme_max;

    }
    if(e[pos].type == right){ //mise a jour des vitesses avec mur "froid"
      int particule = e[pos].ia;

      double norme = sqrt(pow(p[particule].vx,2)+pow(p[particule].vy,2));
      double rd = drand48();

      //on norme la vitesse
      p[particule].vx = p[particule].vx/norme;
      p[particule].vy = p[particule].vy/norme;

      //on modelise maxwell-boltzmann par une loi de proba uniforme
      p[particule].vx *= -rd*norme_min;
      p[particule].vy *= rd*norme_min;
    }


    if(e[pos].type == top||e[pos].type == bottom){
      int particule = e[pos].ia;
      p[particule].vy *= -1;
    }
    
    //code pour rebond elastique avec les murs
    /*
    if(e[pos].type == left||e[pos].type == right){
      int particule = e[pos].ia;
      p[particule].vx *= -1;
    }
    */
    
}


void wall_bounce(Particle *p, int Np, double Lmin, double Lmax, double diameter, double lgmax){
  int i;

  //confinement de la paroi
  if(p[Np].x > lgmax-pow(10,-14)){
    p[Np].x = lgmax-pow(10,-14);
    p[Np].vx *= -1;
  }
  if(p[Np].x < Lmin+pow(10,-14)){
    p[Np].x = Lmin+pow(10,-14);
    p[Np].vx *= -1;
  }

  //confinement des particules
  for(i=0; i<Np; i++){
    if(p[i].x <= Lmin){
      p[i].x = Lmin +diameter/2 +pow(10,-14);
    }
    if(p[i].x >= lgmax){
      p[i].x = lgmax - diameter/2 -pow(10,-14);
    }
    if(p[i].y <=Lmin){
      p[i].y = Lmin +diameter/2 +pow(10,-14);
    }
    if(p[i].y >= Lmax){
      p[i].y = Lmax -diameter/2 -pow(10,-14);
    }
    if(p[i].cote ==0 && p[i].x>=p[Np].x){ //gauche
      p[i].x = p[Np].x -diameter/2 -pow(10,-14);
    }
    if(p[i].cote ==1 && p[i].x<=p[Np].x){ //droite
      p[i].x = p[Np].x +diameter/2 +pow(10,-14);
    }
  }
}


void initparticles( Particle *p, int Np, double Lmin, double Lmax, double diameter, double lgmax){
  int i , d;
  double dx , dy;

  for( i=0;i<Np;i++){ 
    p[i].cote = rand()%2; //chaque particule est soit à gauche (cote=0) soit soit à droite(cote=1)
    

    if(p[i].cote == 0){ //si particule est gauche
      p[i].x = Lmin +diameter/2 + (lgmax/2-Lmin-diameter)*drand48(); //random positions for intial condition
      p[i].y = Lmin +diameter/2 + (Lmax-Lmin-diameter)*drand48();
      
      //placement des particules pour eviter overlap 
      if(i>0){ 
        for(int j=0; j<i ; j++){
          dx = p[i].x-p[j].x;
          dy = p[i].y-p[j].y;
          d = sqrt(dx*dx + dy*dy);

          if(d<diameter){
            p[i].x = Lmin +diameter/2 + (lgmax/2-Lmin-diameter)*drand48(); //random positions for intial condition
            p[i].y = Lmin +diameter/2 + (Lmax-Lmin-diameter)*drand48();
            j=-1; //recommencer si particule mal placée
          }
        }
      }
      
      
    }
    if(p[i].cote == 1){ //si particule est droite
      p[i].x = lgmax/2 +diameter/2 + (lgmax-lgmax/2-diameter)*drand48(); //random positions for intial condition
      p[i].y = Lmin +diameter/2 + (Lmax-Lmin-diameter)*drand48();

      //placement des particules pour eviter overlap
      if(i>0){
        for(int j=0; j<i ; j++){
          dx = p[i].x-p[j].x;
          dy = p[i].y-p[j].y;
          d = sqrt(dx*dx + dy*dy);

          if(d<diameter){
            p[i].x = lgmax/2 +diameter/2 + (lgmax-lgmax/2-diameter)*drand48(); //random positions for intial condition
            p[i].y = Lmin +diameter/2 + (Lmax-Lmin-diameter)*drand48();
            j=-1;//recommencer si particule mal placée
          }
        }
      }
    }
    p[i].vx = 2*(drand48()-0.5);// choose random speeds too using drand48();
    p[i].vy = 2*(drand48()-0.5);
  }

  //position initiale de la paroi
  p[Np].x = lgmax/2; //on place la paroi au milieu de l'enceinte
  p[Np].y = Lmin;
  p[Np].vx = 0;
  p[Np].vy = 0;
}

int main(){
  double T_ech = 0.125;  //T_ech =1/8 pour eviter problemes de derive
  double FPS=50; 
  int Np=200; //nombre de particules (le programme ne tournait pas au dessus de 200 particules)
  double diameter=0.2;
  int Pix=450; 
  double Lmax=10, Lmin=0; 
  double lgmax = 3*Lmax; //grande longeur de la boite
  double tmin;
  int Tmax = 20000; //durée du programme

  Mem_paroi M;

  double f = 1; //coeff de frottement f=1 pas de frottements

  //tracé histogramme
  int nb_chocs = 0;
  double vitesse_max = 3; //arbitraire, 3 est une bonne valeur ici
  int Taille_Hist = 200; //tranches de l'histogramme
  int Hist[Taille_Hist];

  Graphics gw(Np,Pix, Lmin ,Lmax,diameter);
  srand48(time(NULL));//inititalize random numbers -- to find always the same value // you can replace "1" by time(NULL) 

  Particle *p= (Particle *) malloc( (Np+1) *sizeof( Particle)); //an array of particles +1 pour la parroi
  initparticles(p,Np,Lmin, Lmax,diameter, lgmax); //place particles in box
  Event *e = (Event *) malloc((4*Np+1+(Np*(Np-1))/2+Np)* sizeof(Event)); //malloc 4Np+1 with wall, Np(Np-1)/2 with particles, Np with paroi
  
  gw.draw(p,FPS,0); //draw initial position and pause one second
  sleep(1);

  //set initial paramter to animations
  e[4*Np].type = animation;
  e[4*Np].ia = -1;
  e[4*Np].time = T_ech;
  
  //set initial collision to marker to particle
  M.type=particle;
  M.ia=-1; 
  M.ib=-1;
  
  for (int l=0; l < Tmax;l++){
    
    //calculs des temps de collision
    compute_parroi_collision_time(p,Np,e,diameter);
    compute_particle_collision_time(p,Np,e,diameter);
    compute_wall_collision_time(p, Np,e, diameter, Lmax, Lmin, lgmax);
    int pos  = compute_min_time(Np,e,M);
    tmin=e[pos].time;

    if(e[pos].type != animation){
      new_pos(p, Np, tmin);
      nb_chocs += 1; //décompte nombre de chocs pour répartition des vitesses
      new_speed(p,pos,e, Np);
      wall_bounce(p,Np,Lmin,Lmax,diameter, lgmax);
      e[4*Np].time -= tmin;
      
      //mise a jour du marker pour eviter les problemes de precision et donc chocs repetitifs avec paroi
      M.type=e[pos].type;
      M.ia=e[pos].ia;
      M.ib=e[pos].ib;
      
    }
    else
    {
      //animation
      new_pos(p, Np, tmin);
      wall_bounce(p,Np,Lmin,Lmax,diameter, lgmax);
      gw.draw(p, FPS,l);
      e[4*Np].time = T_ech;

      //double p_x = new_pos(p,Np,tmin); //decommenter si on veut recuperer abscisse de la paroi en fonction du temps
      //printf("%d %f\n",l,p_x);
    } 

    if(nb_chocs>100){ //remplissage histogramme
      //printf("%d\n", nb_chocs);
      for(int i=0;i<Np;i++){
        double norme = sqrt(p[i].vx*p[i].vx + p[i].vy*p[i].vy)/vitesse_max;
        if(norme*Taille_Hist < Taille_Hist){
          Hist[(int)(norme*Taille_Hist)] +=1;
        }
        /*else{
            printf("%f>%d\n", norme*Taille_Hist, Taille_Hist);
        }*/
      }
    }//fin remplissage Hist
  }//fin boucle for

  
  //ouverture du fichier
  FILE*file2=NULL;
  file2=fopen("Distribution.txt","w"); //Nom du fichier
  if(file2 == NULL){
    printf("Error: cannot open file");
    exit(1);
  }


  //historgramme centré
  double widthBins = vitesse_max/Taille_Hist;
  double widthBinsHalf = widthBins/2;
  for(int i=0; i<Taille_Hist; i++){
    fprintf(file2, "%f %d\n", i*widthBins+widthBinsHalf, Hist[i]); //recentrer les normes
  }

  fflush(file2);
  pclose(file2);

  free(p);
  free(e);

  return 0;
}
