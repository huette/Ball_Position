// Ball_Position.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//

#include "stdafx.h"
#include <math.h>
#include <stdio.h>


typedef struct vec2d vec2d;

struct vec2d {
        float x;
        float y;
        
};


typedef struct player  player;

struct player {
        vec2d  p;
        float r;
        
};
 
/* Return the difference of two vectors, (vector1 - vector2). */
float vdistance(const vec2d vector1, const vec2d vector2)
{
       float d;
      d = sqrtf(powf(vector1.x - vector2.x,2) + powf(vector1.y - vector2.y, 2));
              
        return d;
}

/* Return the difference of two vectors, (vector1 - vector2). */
vec2d vdiff(const vec2d vector1, const vec2d vector2)
{
        vec2d v;
        v.x = vector1.x - vector2.x;
        v.y = vector1.y - vector2.y;
      
        return v;
}

/* Return the difference of two vectors, (vector1 - vector2). */
vec2d vdiffmean(const vec2d vector1, const vec2d vector2)
{
        vec2d v;
        v.x = (vector1.x + vector2.x)/2;
        v.y = (vector1.y + vector2.y)/2;
      
        return v;
}
 
/* Return the sum of two vectors. */
vec2d vsum(const vec2d vector1, const vec2d vector2)
{
        vec2d v;
        v.x = vector1.x + vector2.x;
        v.y = vector1.y + vector2.y;
      
        return v;
}
 
/* Multiply vector by a number. */
vec2d vmul(const vec2d vector, const double n)
{
        vec2d v;
        v.x = vector.x * n;
        v.y = vector.y * n;
       
        return v;
}
 
/* Divide vector by a number. */
vec2d vdiv(const vec2d vector, const double n)
{
        vec2d v;
        v.x = vector.x / n;
        v.y = vector.y / n;
      
        return v;
}
 
/* Return the Euclidean norm. */
double vnorm(const vec2d vector)
{
        return sqrt(vector.x * vector.x + vector.y * vector.y );
}
 
/* Return the dot product of two vectors. */
double dot(const vec2d vector1, const vec2d vector2)
{
        return vector1.x * vector2.x + vector1.y * vector2.y ;
}
 
/* Replace vector with its cross product with another vector. 
vec2d cross(const vec2d vector1, const vec2d vector2)
{
        vec2d v;
        v.x = vector1.y * vector2.z - vector1.z * vector2.y;
        v.y = vector1.z * vector2.x - vector1.x * vector2.z;
        v.z = vector1.x * vector2.y - vector1.y * vector2.x;
        return v;
}
*/
double cross(const vec2d vector1, const vec2d vector2)
{
       
        return vector1.x * vector2.y - vector1.y * vector2.x;
       
        
}
int intersection(vec2d *result1, vec2d *result2,
                  const vec2d p1, const float r1,
                  const vec2d p2, const float r2)

{
  float a, dx, dy, d, h, rx, ry;
  float x2, y2;
  vec2d t1, t2;

  /* dx and dy are the vertical and horizontal distances between
   * the circle centers.
   */
  dx = p2.x- p1.x;
  dy = p2.y - p1.y;

  /* Determine the straight-line distance between the centers. */
  d = sqrt((dy*dy) + (dx*dx));
  //d = hypot(dx,dy); // Suggested by Keith Briggs

  /* Check for solvability. */
  if (d > (r1 + r2))
  {
    /* no solution. circles do not intersect. */
	  t1.x = dx;
	  t1.y = dy;
	  t1 = vmul(t1, r1/d);
	  t2 = vsum(p1, t1);
	  *result1 = t1;
	  *result2 = t2;
    return 0;
  }
  if (d < fabs(r1 - r2))
  {
    /* no solution. one circle is contained in the other */
    return 0;
  }

  /* 'point 2' is the point where the line through the circle
   * intersection points crosses the line between the circle
   * centers.  
   */

  /* Determine the distance from point 0 to point 2. */
  a = ((r1*r1) - (r2*r2) + (d*d)) / (2.0 * d) ;

  /* Determine the coordinates of point 2. */
  x2 = p1.x + (dx * a/d);
  y2 = p1.y + (dy * a/d);

  /* Determine the distance from point 2 to either of the
   * intersection points.
   */
  h = sqrt((r1*r1) - (a*a));

  /* Now determine the offsets of the intersection points from
   * point 2.
   */
  rx = -dy * (h/d);
  ry = dx * (h/d);

  /* Determine the absolute intersection points. */
   t1.x = x2 + rx;
  t2.x = x2 - rx;
  t1.y = y2 + ry;
  t2.y = y2 - ry;

  *result1=t1;
  *result2=t2;
  

  return 1;
}

float triangleArea(const vec2d p1,const vec2d p2, const vec2d p3) {
	
	float a,b,c,s, area;
	a = vdistance(p1,p2);
	b = vdistance(p2,p3);
	c = vdistance(p3,p1);
	s = (a+b +c)/2;
	area = sqrtf(s*(s-a)*(s-b)*(s-c));
	return area;
}

int _tmain(int argc, _TCHAR* argv[])
{
	player players[10];
	vec2d solutions[5];
     float dis1, dis2, dis3, dis4, nearest, area, maxarea;
      int     result, min, max, numPlayer, thirdpoint;
	
	numPlayer = 5;
	min = 0;
	max = 0;
	maxarea = 0;

	players[0].p.x = 2;
	players[0].p.y = 2;
	players[0].r = 1.41;

	players[1].p.x = 4;
	players[1].p.y = 4;
	players[1].r = 1.55;

	players[2].p.x = 4.5;
	players[2].p.y = 2;
	players[2].r = 1.78;

	players[3].p.x = 2.4;
	players[3].p.y = 2;
	players[3].r = 1.16;

	players[4].p.x = 3.5;
	players[4].p.y = 2.5;
	players[4].r = 0.7;

	if (numPlayer==2) {
		result = intersection(&solutions[0], &solutions[1], players[0].p, players[0].r, players[1].p, players[1].r);
		solutions[4] = vdiffmean(solutions[0], solutions[1]);
	}
	else 
	{

	
	
	for (int i=0; i<numPlayer; i++) {
		
		if (players[i].r<players[min].r) min=i;
		if (players[i].r>players[max].r) max=i;
	}
	
	
	for (int i=0; i<numPlayer;i++) {
		if (i!=min && i!=max) {
			area = triangleArea(players[min].p,players[max].p, players[i].p);
			if (area>maxarea) 
				{ 
					thirdpoint = i;
					maxarea = area;
				}
		}
	
	}
	printf("Triangle %d - %d - %d \n", min, max, thirdpoint);

				result = intersection(&solutions[0], &solutions[1], players[min].p, players[min].r, players[max].p, players[max].r);
		
				result = intersection(&solutions[2], &solutions[3], players[min].p, players[min].r, players[thirdpoint].p, players[thirdpoint].r);
		
				dis1 = vdistance(solutions[0], solutions[2]);
				nearest = dis1;
				solutions[4] = vdiffmean(solutions[0], solutions[2]);
				
				dis2 = vdistance(solutions[0], solutions[3]);
				dis3 = vdistance(solutions[1], solutions[2]);
				dis4 = vdistance(solutions[1], solutions[3]);

				if (dis2<nearest) {
						nearest = dis2;
						solutions[4] = vdiffmean(solutions[0], solutions[3]);
				}
				if (dis3<nearest) {
						nearest = dis3;
						solutions[4] = vdiffmean(solutions[1], solutions[2]);
				}
				if (dis4<nearest) {
						nearest = dis4;
						solutions[4] = vdiffmean(solutions[1], solutions[3]);
				}
	}	
	printf("Solution x:%g y: %g \n",solutions[4].x, solutions[4].y);
}

