// Mars lander simulator
// Version 1.8
// Graphics functions
// Gabor Csanyi and Andrew Gee, October 2014

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

// Some functions adapted from freeglut_geometry.c, which is covered by the
// following license:
//
// Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
// Written by Pawel W. Olszta, <olszta@sourceforge.net>
// Creation date: Fri Dec 3 1999
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
// PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
// IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// Some functions adapted from trackball.cpp by Gavin Bell, which is covered by
// the following license:
//
// (c) Copyright 1993, 1994, Silicon Graphics, Inc.
// ALL RIGHTS RESERVED
// Permission to use, copy, modify, and distribute this software for
// any purpose and without fee is hereby granted, provided that the above
// copyright notice appear in all copies and that both the copyright notice
// and this permission notice appear in supporting documentation, and that
// the name of Silicon Graphics, Inc. not be used in advertising
// or publicity pertaining to distribution of the software without specific,
// written prior permission.
//
// THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
// AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
// INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
// FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
// GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
// SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
// KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
// LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
// THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
// ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
// POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
//
// US Government Users Restricted Rights
// Use, duplication, or disclosure by the Government is subject to
// restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
// (c)(1)(ii) of the Rights in Technical Data and Computer Software
// clause at DFARS 252.227-7013 and/or in similar or successor
// clauses in the FAR or the DOD or NASA FAR Supplement.
// Unpublished-- rights reserved under the copyright laws of the
// United States.  Contractor/manufacturer is Silicon Graphics,
// Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
//
// OpenGL(TM) is a trademark of Silicon Graphics, Inc.

#define DECLARE_GLOBAL_VARIABLES
#include "lander.h"

void invert (double m[], double mout[])
  // Inverts a 4x4 OpenGL rotation matrix
{
  double zero_three, one_three, two_three;
  zero_three = -m[12]*m[0] - m[13]*m[1] - m[14]*m[2];
  one_three = -m[12]*m[4] - m[13]*m[5] - m[14]*m[6];
  two_three = -m[12]*m[8] - m[13]*m[9] - m[14]*m[10];
  mout[1] = m[4]; mout[4] = m[1]; mout[2] = m[8]; mout[8] = m[2]; 
  mout[6] = m[9]; mout[9] = m[6]; mout[12] = zero_three; mout[13] = one_three; 
  mout[14] = two_three; mout[0] = m[0]; mout[5] = m[5]; mout[10] = m[10];
  mout[15] = 1.0; mout[3] = 0.0; mout[7] = 0.0; mout[11] = 0.0;
}

void xyz_euler_to_matrix (vector3d ang, double m[])
  // Constructs a 4x4 OpenGL rotation matrix from xyz Euler angles
{
  double sin_a, sin_b, sin_g, cos_a, cos_b, cos_g;
  double ra, rb, rg;

  // Pre-calculate radian angles
  ra = ang.x*M_PI/(double)180;
  rb = ang.y*M_PI/(double)180;
  rg = ang.z*M_PI/(double)180;

  // Pre-calculate sines and cosines
  cos_a = cos(ra);
  cos_b = cos(rb);
  cos_g = cos(rg);
  sin_a = sin(ra);
  sin_b = sin(rb);
  sin_g = sin(rg);

  // Create the correct matrix coefficients
  m[0] = cos_a * cos_b;
  m[1] = sin_a * cos_b;
  m[2] = - sin_b;
  m[3] = 0.0;
  m[4] = cos_a * sin_b * sin_g - sin_a * cos_g;
  m[5] = sin_a * sin_b * sin_g + cos_a * cos_g;
  m[6] = cos_b * sin_g;
  m[7] = 0.0;
  m[8] = cos_a * sin_b * cos_g + sin_a * sin_g;
  m[9] = sin_a * sin_b * cos_g - cos_a * sin_g;
  m[10] = cos_b * cos_g;
  m[11] = 0.0;
  m[12] = 0.0;
  m[13] = 0.0;
  m[14] = 0.0;
  m[15] = 1.0;
}

vector3d matrix_to_xyz_euler (double m[])
  // Decomposes a 4x4 OpenGL rotation matrix into xyz Euler angles
{
  double tmp;
  vector3d ang;

  // Catch degenerate elevation cases
  if (m[2] < -0.99999999) {
    ang.y = 90.0;
    ang.x = 0.0;
    ang.z = acos(m[8]);
    if ( (sin(ang.z)>0.0) ^ (m[4]>0.0) ) ang.z = -ang.z;
    ang.z *= 180.0/M_PI;
    return ang;
  }
  if (m[2] > 0.99999999) {
    ang.y = -90.0;
    ang.x = 0.0;
    ang.z = acos(m[5]);
    if ( (sin(ang.z)<0.0) ^ (m[4]>0.0) ) ang.z = -ang.z;
    ang.z *= 180.0/M_PI;
    return ang;
  }

  // Non-degenerate elevation - between -90 and +90
  ang.y = asin( -m[2] );

  // Now work out azimuth - between -180 and +180
  tmp = m[0]/cos(ang.y); // the denominator will not be zero
  if ( tmp <= -1.0 ) ang.x = M_PI;
  else if ( tmp >= 1.0 ) ang.x = 0.0;
  else ang.x = acos( tmp );
  if ( ((sin(ang.x) * cos(ang.y))>0.0) ^ ((m[1])>0.0) ) ang.x = -ang.x;

  // Now work out roll - between -180 and +180
  tmp = m[10]/cos(ang.y); // the denominator will not be zero
  if ( tmp <= -1.0 ) ang.z = M_PI;
  else if ( tmp >= 1.0 ) ang.z = 0.0;
  else ang.z = acos( tmp );
  if ( ((sin(ang.z) * cos(ang.y))>0.0) ^ ((m[6])>0.0) ) ang.z = -ang.z;

  // Convert to degrees
  ang.y *= 180.0/M_PI;
  ang.x *= 180.0/M_PI;
  ang.z *= 180.0/M_PI;

  return ang;
}

void normalize_quat (quat_t &q)
  // Normalizes a quaternion
{
  double mag;
  mag = (q.v.x*q.v.x + q.v.y*q.v.y + q.v.z*q.v.z + q.s*q.s);
  if (mag > 0.0) {
    q.v.x /= mag; q.v.y /= mag;
    q.v.z /= mag; q.s /= mag;
  }
}

quat_t axis_to_quat (vector3d a, const double phi)
  // Given an axis and angle, compute quaternion
{
  quat_t q;
  q.v = a.norm() * sin(phi/2.0);
  q.s = cos(phi/2.0);
  return q;
}

double project_to_sphere (const double r, const double x, const double y)
  // Project an x,y pair onto a sphere of radius r or a hyperbolic sheet if
  // we are away from the centre of the sphere
{
  double d, t, z;

  d = sqrt(x*x + y*y);
  if (d < (r * 0.70710678118654752440)) z = sqrt(r*r - d*d);
  else { // on hyperbola
    t = r / 1.41421356237309504880;
    z = t*t / d;
  }
  return z;
}

quat_t add_quats (quat_t q1, quat_t q2)
  // Given two rotations q1 and q2, calculate single equivalent quaternion
{
  quat_t s;
  vector3d t1 = q1.v * q2.s;
  vector3d t2 = q2.v * q1.s;
  vector3d t3 = q2.v ^ q1.v;
  vector3d tf = t1+t2;
  s.v = tf+t3;
  s.s = q1.s * q2.s - q1.v * q2.v;
  normalize_quat(s);
  return s;
}

void quat_to_matrix (double m[], const quat_t q)
  // Convert quaternion into a rotation matrix
{
  m[0] = 1.0 - 2.0 * (q.v.y * q.v.y + q.v.z * q.v.z);
  m[1] = 2.0 * (q.v.x * q.v.y - q.v.z * q.s);
  m[2] = 2.0 * (q.v.z * q.v.x + q.v.y * q.s);
  m[3] = 0.0;

  m[4] = 2.0 * (q.v.x * q.v.y + q.v.z * q.s);
  m[5] = 1.0 - 2.0 * (q.v.z * q.v.z + q.v.x * q.v.x);
  m[6] = 2.0 * (q.v.y * q.v.z - q.v.x * q.s);
  m[7] = 0.0;

  m[8] = 2.0 * (q.v.z * q.v.x - q.v.y * q.s);
  m[9] = 2.0 * (q.v.y * q.v.z + q.v.x * q.s);
  m[10] = 1.0 - 2.0 * (q.v.y * q.v.y + q.v.x * q.v.x);
  m[11] = 0.0;

  m[12] = 0.0; m[13] = 0.0; m[14] = 0.0; m[15] = 1.0;
}

quat_t track_quats (const double p1x, const double p1y, const double p2x, const double p2y)
  // Derive quaternion from x and y mouse displacements
{
  double t, phi;
  vector3d a, p1, p2, d;
  quat_t q;

  if ((p1x == p2x) && (p1y == p2y)) {
    q.v.x = 0.0; q.v.y = 0.0; q.v.z = 0.0; q.s = 1.0;
    return q;
  }

  p1.x = p1x; p1.y = p1y;
  p1.z = project_to_sphere(0.5, p1x, p1y);
  p2.x = p2x; p2.y = p2y;
  p2.z = project_to_sphere(0.5, p2x, p2y);
  a = p2^p1; d = p1-p2; t = d.abs();
  if (t > 1.0) t = 1.0;
  if (t < -1.0) t = -1.0;
  phi = 2.0 * asin(t);
  return axis_to_quat(a, phi);
}

void microsecond_time (unsigned long long &t)
  // Returns system time in microseconds, used for introducing delays at low simulation speeds
{
#ifdef WIN32
  LARGE_INTEGER counter, frequency;
  QueryPerformanceFrequency(&frequency);
  frequency.QuadPart /= 1000000;
  QueryPerformanceCounter(&counter);
  t = (unsigned long long)(counter.QuadPart/frequency.QuadPart);
#else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  t = (unsigned long long)tv.tv_usec + 1000000 * (unsigned long long)tv.tv_sec;
#endif
}

void fghCircleTable (double **sint, double **cost, const int n)
  // Borrowed from freeglut source code, used to draw hemispheres and open cones
{
  int i;
  const int size = abs(n);
  const double angle = 2*M_PI/(double)( ( n == 0 ) ? 1 : n );
  
  *sint = (double*) calloc(sizeof(double), size+1);
  *cost = (double*) calloc(sizeof(double), size+1);
  if (!(*sint) || !(*cost)) exit(1);
  
  (*sint)[0] = 0.0;
  (*cost)[0] = 1.0;
  
  for (i=1; i<size; i++) {
    (*sint)[i] = sin(angle*i);
    (*cost)[i] = cos(angle*i);
  }
  
  (*sint)[size] = (*sint)[0];
  (*cost)[size] = (*cost)[0];
}

void glutOpenHemisphere (GLdouble radius, GLint slices, GLint stacks)
  // Modified from freeglut's glutSolidSphere
{
  int i, j;
  double z0, z1, r0, r1, *sint1, *cost1, *sint2, *cost2;
  
  fghCircleTable(&sint1, &cost1, -slices);
  fghCircleTable(&sint2, &cost2, stacks*2);
  z1 = cost2[(stacks>0)?1:0];
  r1 = sint2[(stacks>0)?1:0];
  
  // Middle stacks
  for (i=1; i<stacks-1; i++) {
    z0 = z1; z1 = cost2[i+1];
    r0 = r1; r1 = sint2[i+1];
    if ((z1 > 0) || (z0 > 0)) continue; // hemisphere
    glBegin(GL_QUAD_STRIP);
    for (j=0; j<=slices; j++) {
      glNormal3d(cost1[j]*r1, sint1[j]*r1, z1);
      glVertex3d(cost1[j]*r1*radius, sint1[j]*r1*radius, z1*radius);
      glNormal3d(cost1[j]*r0, sint1[j]*r0, z0);
      glVertex3d(cost1[j]*r0*radius, sint1[j]*r0*radius, z0*radius);
    }
    glEnd();
  }

  // Bottom cap
  z0 = z1; r0 = r1;
  glBegin(GL_TRIANGLE_FAN);
  glNormal3d(0,0,-1);
  glVertex3d(0,0,-radius);
  for (j=0; j<=slices; j++) {
    glNormal3d(cost1[j]*r0, sint1[j]*r0, z0);
    glVertex3d(cost1[j]*r0*radius, sint1[j]*r0*radius, z0*radius);
  }
  glEnd();
  
  free(sint1); free(cost1);
  free(sint2); free(cost2);
}

void glutMottledSphere (GLdouble radius, GLint slices, GLint stacks)
  // Modified from freeglut's glutSolidSphere, we use this to draw a mottled sphere by modulating
  // the vertex colours.
{
    int i, j;
    unsigned short rtmp = 0;
    double z0, z1, r0, r1, *sint1, *cost1, *sint2, *cost2;
    double *rnd1, *rnd2, *new_r, *old_r, *tmp;
    double mottle = 0.2;

    fghCircleTable(&sint1, &cost1, -slices);
    fghCircleTable(&sint2, &cost2, stacks*2);
    rnd1 = (double*) calloc(sizeof(double), slices+1);
    rnd2 = (double*) calloc(sizeof(double), slices+1);
    z0 = 1.0; z1 = cost2[(stacks>0)?1:0];
    r0 = 0.0; r1 = sint2[(stacks>0)?1:0];

    // Top cap
    glBegin(GL_TRIANGLE_FAN);
    glNormal3d(0,0,1);
    glColor3f(0.63, 0.33, 0.22);
    glVertex3d(0,0,radius);
    new_r = rnd1;
    for (j=slices; j>=0; j--) {
      glNormal3d(cost1[j]*r1, sint1[j]*r1, z1);
      if (j) {
	new_r[j] = (1.0-mottle) + mottle*randtab[rtmp];
	rtmp = (rtmp+1)%N_RAND;
      } else new_r[j] = new_r[slices];
      glColor3f(new_r[j]*0.63, new_r[j]*0.33, new_r[j]*0.22);
      glVertex3d(cost1[j]*r1*radius, sint1[j]*r1*radius, z1*radius);
    }
    glEnd();

    // Middle stacks
    old_r = rnd1; new_r = rnd2;
    for (i=1; i<stacks-1; i++) {
      z0 = z1; z1 = cost2[i+1];
      r0 = r1; r1 = sint2[i+1];
      glBegin(GL_QUAD_STRIP);
      for (j=0; j<=slices; j++) {
	glNormal3d(cost1[j]*r1, sint1[j]*r1, z1);
	if (j != slices) {
	  new_r[j] = (1.0-mottle) + mottle*randtab[rtmp];
	  rtmp = (rtmp+1)%N_RAND;
	} else new_r[j] = new_r[0];
	glColor3f(new_r[j]*0.63, new_r[j]*0.33, new_r[j]*0.22);
	glVertex3d(cost1[j]*r1*radius, sint1[j]*r1*radius, z1*radius);
	glNormal3d(cost1[j]*r0, sint1[j]*r0, z0);
	glColor3f(old_r[j]*0.63, old_r[j]*0.33, old_r[j]*0.22);
	glVertex3d(cost1[j]*r0*radius, sint1[j]*r0*radius, z0*radius);
      }
      tmp = old_r; old_r = new_r; new_r = tmp;
      glEnd();
    }

    // Bottom cap
    z0 = z1; r0 = r1;
    glBegin(GL_TRIANGLE_FAN);
    glNormal3d(0,0,-1);
    glColor3f(0.63, 0.33, 0.22);
    glVertex3d(0,0,-radius);
    for (j=0; j<=slices; j++) {
      glNormal3d(cost1[j]*r0, sint1[j]*r0, z0);
      glColor3f(old_r[j]*0.63, old_r[j]*0.33, old_r[j]*0.22);
      glVertex3d(cost1[j]*r0*radius, sint1[j]*r0*radius, z0*radius);
    }
    glEnd();

    free(rnd1); free(rnd2);
    free(sint1); free(cost1);
    free(sint2); free(cost2);
}

void glutCone (GLdouble base, GLdouble height, GLint slices, GLint stacks, bool closed)
  // Modified from freeglut's glutSolidCone, we need this (a) to draw cones without bases and
  // (b) to draw cones with bases, which glutSolidCone does not do correctly under Windows,
  // for some reason.
{
  int i, j;
  double z0, z1, r0, r1, *sint, *cost;
  const double zStep = height / ( ( stacks > 0 ) ? stacks : 1 );
  const double rStep = base / ( ( stacks > 0 ) ? stacks : 1 );
  const double cosn = ( height / sqrt ( height * height + base * base ));
  const double sinn = ( base   / sqrt ( height * height + base * base ));
  
  fghCircleTable(&sint, &cost, -slices);
  z0 = 0.0; z1 = zStep;
  r0 = base; r1 = r0 - rStep;
  
  if (closed) {
    glBegin(GL_TRIANGLE_FAN);
    glNormal3d(0.0, 0.0, -1.0);
    glVertex3d(0.0, 0.0, z0);
    for (j=0; j<=slices; j++) glVertex3d(cost[j]*r0, sint[j]*r0, z0);
    glEnd();
  }

  for (i=0; i<stacks-1; i++) {
    glBegin(GL_QUAD_STRIP);
    for (j=0; j<=slices; j++) {
      glNormal3d(cost[j]*sinn, sint[j]*sinn, cosn);
      glVertex3d(cost[j]*r0, sint[j]*r0, z0);
      glVertex3d(cost[j]*r1, sint[j]*r1, z1);
    }
    z0 = z1; z1 += zStep;
    r0 = r1; r1 -= rStep;
    glEnd();
  }
  
  glBegin(GL_TRIANGLES);
  glNormal3d(cost[0]*sinn, sint[0]*sinn, cosn);
  for (j=0; j<slices; j++) {
    glVertex3d(cost[j+0]*r0, sint[j+0]*r0, z0);
    glVertex3d(0.0, 0.0, height);
    glNormal3d(cost[j+1]*sinn, sint[j+1]*sinn, cosn);
    glVertex3d(cost[j+1]*r0, sint[j+1]*r0, z0);
  }
  glEnd();
  
  free(sint); free(cost);
}

void enable_lights (void) 
  // Enable the appropriate subset of lights
{
  if (static_lighting) {
    glDisable(GL_LIGHT0); glDisable(GL_LIGHT1);
    glEnable(GL_LIGHT2); glEnable(GL_LIGHT3);
    glDisable(GL_LIGHT4); glDisable(GL_LIGHT5);
  } else {
    glEnable(GL_LIGHT0); glEnable(GL_LIGHT1);
    glDisable(GL_LIGHT2); glDisable(GL_LIGHT3);
    glDisable(GL_LIGHT4); glDisable(GL_LIGHT5);
  }
}

void setup_lights (void)
  // Specifies attributes of all lights, enables a subset of lights according to the lighting model
{
  GLfloat none[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat low[] = { 0.15, 0.15, 0.15, 1.0 };
  GLfloat medium[] = { 0.5, 0.5, 0.5, 1.0 };
  GLfloat high[] = { 0.75, 0.75, 0.75, 1.0 };

  // Lights 0 and 1 are for the dynamic lighting model, with the lights fixed in the viewer's reference frame
  glLightfv(GL_LIGHT0, GL_AMBIENT, none);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, high);
  glLightfv(GL_LIGHT0, GL_SPECULAR, none);
  glLightfv(GL_LIGHT0, GL_POSITION, top_right);
  glLightfv(GL_LIGHT1, GL_AMBIENT, none);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, medium);
  glLightfv(GL_LIGHT1, GL_SPECULAR, none);
  glLightfv(GL_LIGHT1, GL_POSITION, straight_on);

  // Lights 2 and 3 are for the static lighting model, with the lights fixed in the planetary reference frame
  glLightfv(GL_LIGHT2, GL_AMBIENT, none);
  glLightfv(GL_LIGHT2, GL_DIFFUSE, high);
  glLightfv(GL_LIGHT2, GL_SPECULAR, none);
  glLightfv(GL_LIGHT3, GL_AMBIENT, low);
  glLightfv(GL_LIGHT3, GL_DIFFUSE, none);
  glLightfv(GL_LIGHT3, GL_SPECULAR, none);

  // Lights 4 and 5 are for highlighting the lander with static lights, to avoid flat views with ambient illumination only
  glLightfv(GL_LIGHT4, GL_AMBIENT, none);
  glLightfv(GL_LIGHT4, GL_DIFFUSE, low);
  glLightfv(GL_LIGHT4, GL_SPECULAR, none);
  glLightfv(GL_LIGHT5, GL_AMBIENT, none);
  glLightfv(GL_LIGHT5, GL_DIFFUSE, low);
  glLightfv(GL_LIGHT5, GL_SPECULAR, none);

  enable_lights();
}

void glut_print (float x, float y, string s)
  // Prints string at location (x,y) in a bitmap font
{
  unsigned short i;

  glRasterPos2f(x, y);
  for (i = 0; i < s.length(); i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, s[i]);
}

double atmospheric_density (vector3d pos)
  // Simple exponential model between surface and exosphere (around 200km), surface density is approximately 0.017 kg/m^3,
  // scale height is approximately 11km
{
  double alt;

  alt = pos.abs()-MARS_RADIUS;
  if ((alt > EXOSPHERE) || (alt < 0.0)) return 0.0;
  else return (0.017 * exp(-alt/11000.0));
}

void draw_dial (double cx, double cy, double val, string title, string units)
  // Draws a single instrument dial, position (cx, cy), value val, title
{
  int i, e;
  double a;
  ostringstream s;

  // Work out value mantissa and exponent
  if (val <= 0.0) {
    e = 0;
    a = 0.0;
  } else {
    e = 0;
    a = val;
    while (a >= 10.0) {
      e++;
      a /= 10.0;
    }
  }

  // Draw dial ticks
  glColor3f(1.0, 1.0, 1.0);
  glBegin(GL_LINES);
  for (i=30; i<=330; i+=30) {
    glVertex2d(cx - OUTER_DIAL_RADIUS * sin(i*M_PI/180.0), cy - OUTER_DIAL_RADIUS * cos(i*M_PI/180.0));
    glVertex2d(cx - INNER_DIAL_RADIUS * sin(i*M_PI/180.0), cy - INNER_DIAL_RADIUS * cos(i*M_PI/180.0));
  }
  glEnd();

  // Draw dial needle
  glColor3f(0.0, 1.0, 1.0);
  glBegin(GL_LINES);
  glVertex2d(cx, cy);
  glVertex2d(cx - INNER_DIAL_RADIUS * sin((a*30+30)*M_PI/180.0), cy - INNER_DIAL_RADIUS * cos((a*30+30)*M_PI/180.0));
  glEnd();

  // Draw exponent indicator, value and title
  s.precision(1);
  glColor3f(1.0, 1.0, 1.0);
  s.str(""); s << "x 10 ^ " << e << " " << units;
  glut_print(cx+10-3.2*s.str().length(), cy+10, s.str());
  glut_print(cx+10-3.2*title.length(), cy-OUTER_DIAL_RADIUS-15, title);
  s.str(""); s << fixed << val << " " << units;
  glut_print(cx+10-3.2*s.str().length(), cy-OUTER_DIAL_RADIUS-30, s.str());

  // Draw tick labels
  for (i=0; i<=10; i++) {
    switch(i) {
    case 0:
      glut_print(cx - OUTER_DIAL_RADIUS * sin((i*30+30)*M_PI/180.0) - 8, cy - OUTER_DIAL_RADIUS * cos((i*30+30)*M_PI/180.0) - 9, "0");
      break;
    case 1:
      glut_print(cx - OUTER_DIAL_RADIUS * sin((i*30+30)*M_PI/180.0) - 7, cy - OUTER_DIAL_RADIUS * cos((i*30+30)*M_PI/180.0) - 6, "1");
      break;
    case 2:
      glut_print(cx - OUTER_DIAL_RADIUS * sin((i*30+30)*M_PI/180.0) - 9, cy - OUTER_DIAL_RADIUS * cos((i*30+30)*M_PI/180.0) - 4, "2");
      break;
    case 3:
      glut_print(cx - OUTER_DIAL_RADIUS * sin((i*30+30)*M_PI/180.0) - 9, cy - OUTER_DIAL_RADIUS * cos((i*30+30)*M_PI/180.0) - 1, "3");
      break;
    case 4:
      glut_print(cx - OUTER_DIAL_RADIUS * sin((i*30+30)*M_PI/180.0) - 8, cy - OUTER_DIAL_RADIUS * cos((i*30+30)*M_PI/180.0) + 3, "4");
      break;
    case 5:
      glut_print(cx - OUTER_DIAL_RADIUS * sin((i*30+30)*M_PI/180.0) - 3, cy - OUTER_DIAL_RADIUS * cos((i*30+30)*M_PI/180.0) + 4, "5");
      break;
    case 6:
      glut_print(cx - OUTER_DIAL_RADIUS * sin((i*30+30)*M_PI/180.0) + 2, cy - OUTER_DIAL_RADIUS * cos((i*30+30)*M_PI/180.0) + 3, "6");
      break;
    case 7:
      glut_print(cx - OUTER_DIAL_RADIUS * sin((i*30+30)*M_PI/180.0) + 4, cy - OUTER_DIAL_RADIUS * cos((i*30+30)*M_PI/180.0), "7");
      break;
    case 8:
      glut_print(cx - OUTER_DIAL_RADIUS * sin((i*30+30)*M_PI/180.0) + 3, cy - OUTER_DIAL_RADIUS * cos((i*30+30)*M_PI/180.0) - 4, "8");
      break;
    case 9:
      glut_print(cx - OUTER_DIAL_RADIUS * sin((i*30+30)*M_PI/180.0) + 3, cy - OUTER_DIAL_RADIUS * cos((i*30+30)*M_PI/180.0) - 6, "9");
      break;
    case 10:
      glut_print(cx - OUTER_DIAL_RADIUS * sin((i*30+30)*M_PI/180.0) + 3, cy - OUTER_DIAL_RADIUS * cos((i*30+30)*M_PI/180.0) - 8, "10");
      break;
    }
  }
}

void draw_control_bar (double tlx, double tly, double val, double red, double green, double blue, string title)
  // Draws control bar, top left (tlx, tly), val (fraction, range 0-1), colour (red, green, blue), title
{
  glColor3f(red, green, blue);
  glBegin(GL_QUADS);
  glVertex2d(tlx+0.5, tly-19.5);
  glVertex2d(tlx+0.5+239.0*val, tly-19.5);
  glVertex2d(tlx+0.5+239.0*val, tly-0.5);
  glVertex2d(tlx+0.5, tly-0.5);
  glEnd();
  glColor3f(1.0, 1.0, 1.0);
  glBegin(GL_LINE_LOOP);
  glVertex2d(tlx, tly-20.0);
  glVertex2d(tlx+240.0, tly-20.0);
  glVertex2d(tlx+240.0, tly);
  glVertex2d(tlx, tly);
  glEnd();
  glut_print(tlx, tly-40, title);
}

void draw_indicator_lamp (double tcx, double tcy, string off_text, string on_text, bool on)
  // Draws indicator lamp, top centre (tcx, tcy), appropriate text and background colour depending on on/off
{
  if (on) glColor3f(0.5, 0.0, 0.0);
  else glColor3f(0.0, 0.5, 0.0);
  glBegin(GL_QUADS);
  glVertex2d(tcx-74.5, tcy-19.5);
  glVertex2d(tcx+74.5, tcy-19.5);
  glVertex2d(tcx+74.5, tcy-0.5);
  glVertex2d(tcx-74.5, tcy-0.5);
  glEnd();
  glColor3f(1.0, 1.0, 1.0);
  glBegin(GL_LINE_LOOP);
  glVertex2d(tcx-75.0, tcy-20.0);
  glVertex2d(tcx+75.0, tcy-20.0);
  glVertex2d(tcx+75.0, tcy);
  glVertex2d(tcx-75.0, tcy);
  glEnd();
  if (on) glut_print(tcx-70.0, tcy-14.0, on_text);
  else glut_print(tcx-70.0, tcy-14.0, off_text);
}

void draw_instrument_window (void)
  // Draws the instruments
{
  ostringstream s;

  s.precision(1);
  glutSetWindow(instrument_window);
  glClear(GL_COLOR_BUFFER_BIT);

  // Draw altimeter
  draw_dial(view_width + GAP - 400, INSTRUMENT_HEIGHT / 2, altitude, "Altitude", "m");

  // Draw auto-pilot lamp
  draw_indicator_lamp (view_width+GAP-400, INSTRUMENT_HEIGHT-18, "Auto-pilot off", "Auto-pilot on", autopilot_enabled);

  // Draw climb rate meter
  if (climb_speed >= 0.0) draw_dial (view_width+GAP-150, INSTRUMENT_HEIGHT/2, landed ? 0.0 : climb_speed, "Climb rate", "m/s");
  else draw_dial (view_width+GAP-150, INSTRUMENT_HEIGHT/2, landed ? 0.0 : -climb_speed, "Descent rate", "m/s");

  // Draw attitude stabilizer lamp
  draw_indicator_lamp (view_width+GAP-150, INSTRUMENT_HEIGHT-18, "Attitude stabilizer off", "Attitude stabilizer on", stabilized_attitude);

  // Draw ground speed meter
  draw_dial (view_width+GAP+100, INSTRUMENT_HEIGHT/2, landed ? 0.0 : ground_speed, "Ground speed", "m/s");

  // Draw parachute lamp
  switch (parachute_status) {
  case NOT_DEPLOYED:
    draw_indicator_lamp (view_width+GAP+100, INSTRUMENT_HEIGHT-18, "Parachute not deployed", "Do not deploy parachute", !safe_to_deploy_parachute());
    break;
  case DEPLOYED:
    draw_indicator_lamp (view_width+GAP+100, INSTRUMENT_HEIGHT-18, "Parachute deployed", "", false);
    break;
  case LOST:
    draw_indicator_lamp (view_width+GAP+100, INSTRUMENT_HEIGHT-18, "", "Parachute lost", true);
    break;
  }

  // Draw speed bar
  draw_control_bar(view_width+GAP+240, INSTRUMENT_HEIGHT-18, simulation_speed/10.0, 0.0, 0.0, 1.0, "Simulation speed");
  
  // Draw digital clock
  glColor3f(1.0, 1.0, 1.0);
  s.str(""); s << "Time " << fixed << simulation_time << " s";
  glut_print(view_width+GAP+400, INSTRUMENT_HEIGHT-58, s.str());
  if (paused) {
    glColor3f(1.0, 0.0, 0.0);
    glut_print(view_width+GAP+338, INSTRUMENT_HEIGHT-32, "PAUSED");
  }

  // Display coordinates
  glColor3f(1.0, 1.0, 1.0);
  s.str(""); s << "x position " << fixed << position.x << " m";
  glut_print(view_width+GAP+240, INSTRUMENT_HEIGHT-97, s.str());
  s.str(""); s << "velocity " << fixed << velocity_from_positions.x << " m/s";
  glut_print(view_width+GAP+380, INSTRUMENT_HEIGHT-97, s.str());
  s.str(""); s << "y position " << fixed << position.y << " m";
  glut_print(view_width+GAP+240, INSTRUMENT_HEIGHT-117, s.str());
  s.str(""); s << "velocity " << fixed << velocity_from_positions.y << " m/s";
  glut_print(view_width+GAP+380, INSTRUMENT_HEIGHT-117, s.str());
  s.str(""); s << "z position " << fixed << position.z << " m";
  glut_print(view_width+GAP+240, INSTRUMENT_HEIGHT-137, s.str());
  s.str(""); s << "velocity " << fixed << velocity_from_positions.z << " m/s";
  glut_print(view_width+GAP+380, INSTRUMENT_HEIGHT-137, s.str());

  // Draw thrust bar
  s.str(""); s << "Thrust " << fixed << thrust_wrt_world().abs() << " N";
  draw_control_bar(view_width+GAP+240, INSTRUMENT_HEIGHT-170, throttle, 1.0, 0.0, 0.0, s.str());

  // Draw fuel bar
  s.str(""); s << "Fuel " << fixed << fuel*FUEL_CAPACITY << " litres";
  if (fuel > 0.5) draw_control_bar(view_width+GAP+240, INSTRUMENT_HEIGHT-242, fuel, 0.0, 1.0, 0.0, s.str());
  else if (fuel > 0.2) draw_control_bar(view_width+GAP+240, INSTRUMENT_HEIGHT-242, fuel, 1.0, 0.5, 0.0, s.str());
  else draw_control_bar(view_width+GAP+240, INSTRUMENT_HEIGHT-242, fuel, 1.0, 0.0, 0.0, s.str());

  // Display simulation status
  if (landed) glColor3f(1.0, 1.0, 0.0);
  else glColor3f(1.0, 1.0, 1.0);
  s.str(""); s << "Scenario " << scenario;
  if (!landed) s << ": " << scenario_description[scenario];
  glut_print(view_width+GAP-488, 17, s.str());
  if (landed) {
    if (altitude < LANDER_SIZE/2.0) glut_print(80, 17, "Lander is below the surface!");
    else {
      s.str(""); s << "Fuel consumed " << fixed << FUEL_CAPACITY*(1.0-fuel) << " litres";
      glut_print(view_width+GAP-427, 17, s.str());
      s.str(""); s << "Descent rate at touchdown " << fixed << -climb_speed << " m/s";
      glut_print(view_width+GAP-232, 17, s.str());
      s.str(""); s << "Ground speed at touchdown " << fixed << ground_speed << " m/s";
      glut_print(view_width+GAP+16, 17, s.str());
    }
  }

  glutSwapBuffers();
}

void display_help_arrows (void)
  // Displays help arrow in close-up view window
{
  double m[16], p[16], x, y, z, s = -closeup_offset/50.0;
  GLint v[4];
  unsigned short i;
  string ss = "surface";

  glDisable(GL_LIGHTING);
  glColor3f(1.0, 1.0, 1.0);
  glLineWidth(1.0);

  // Surface arrow
  glBegin(GL_LINES);
  glVertex3d(0.0, 2.0*s, 0.0);
  glVertex3d(0.0, 6.0*s, 0.0);
  glEnd();
  glPushMatrix();
  glRotated(-closeup_yr, 0.0, 1.0, 0.0);
  glTranslated(0.0, 6.0*s, 0.0);
  glRotated(90.0, 1.0, 0.0, 0.0);
  glutCone(-0.2*s, -0.5*s, 5, 5, true);
  glGetDoublev(GL_MODELVIEW_MATRIX, m);
  glGetDoublev(GL_PROJECTION_MATRIX, p);
  glGetIntegerv(GL_VIEWPORT, v);
  gluProject(0.0, 0.0, -0.5*s, m, p, v, &x, &y, &z);
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, view_width, 0, view_height, 0.0, 1.0); 
  glRasterPos3f(x-16, y-15, -z);
  for (i = 0; i < ss.length(); i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, ss[i]);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  // Ground speed arrow
  if ((ground_speed > MAX_IMPACT_GROUND_SPEED) && !landed) {
    glBegin(GL_LINES);
    glVertex3d(-2.0*s, 0.0, 0.0);
    glVertex3d(-6.0*s, 0.0, 0.0);
    glEnd();
    glPushMatrix();
    glTranslated(-6.0*s, 0.0, 0.0);
    glRotated(90.0, 0.0, 1.0, 0.0);
    glutCone(-0.2*s, -0.5*s, 5, 5, true);
    glRotated(-90.0, 0.0, 1.0, 0.0);
    glut_print(0.0, 1.25*s, "ground speed");
    glPopMatrix();
  }

  glEnable(GL_LIGHTING);
}

void display_help_text (void)
  // Displays help information in orbital view window
{
  ostringstream s;
  unsigned short i, j;

  glColor3f(1.0, 1.0, 1.0);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, view_width, 0, view_height, -1.0, 1.0); 
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  glut_print(20, view_height-20, "Left arrow - decrease simulation speed");
  glut_print(20, view_height-35, "Right arrow - increase simulation speed");
  glut_print(20, view_height-50, "Space - single step through simulation");

  glut_print(20, view_height-70, "Up arrow - more thrust");
  glut_print(20, view_height-85, "Down arrow - less thrust");

  glut_print(20, view_height-105, "Keys 0-9 - restart simulation in scenario n");

  glut_print(20, view_height-125, "Left mouse - rotate 3D views");
  glut_print(20, view_height-140, "Middle/shift mouse or up wheel - zoom in 3D views");
  glut_print(20, view_height-155, "Right mouse or down wheel - zoom out 3D views");

  glut_print(20, view_height-175, "s - toggle attitude stabilizer");
  glut_print(20, view_height-190, "p - deploy parachute");
  glut_print(20, view_height-205, "a - toggle autopilot");

  glut_print(20, view_height-225, "l - toggle lighting model");
  glut_print(20, view_height-240, "t - toggle terrain texture");
  glut_print(20, view_height-255, "h - toggle help");
  glut_print(20, view_height-270, "Esc/q - quit");

  j = 0;
  for (i=0; i<10; i++) {
    s.str("");
    s << "Scenario " << i << ": " << scenario_description[i];
    if (view_height > 448) glut_print(20, (448-view_height) + view_height-300-15*j, s.str());
    else glut_print(20, view_height-300-15*j, s.str());
    j++;
  }

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

void display_help_prompt (void)
  // Displays help prompt in the close-up view window
{
  string ss = "press 'h' for help";
  unsigned short i;
  unsigned long long t;
  float c;

  microsecond_time(t);
  c = 1.0 - (t - time_program_started) / 3000000.0;
  if (c < 0.0) return;
  glColor4f(1.0, 1.0, 1.0, c);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, view_width, 0, view_height, -1.0, 1.0); 
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);

  glRasterPos2f(view_width/2 - 87, view_height-130);
  for (i = 0; i < ss.length(); i++) glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, ss[i]);

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

void draw_orbital_window (void)
  // Draws the orbital view
{
  unsigned short i, j;
  double m[16], sf;
  GLint slices, stacks;

  glutSetWindow(orbital_window);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Viewing transformation
  quat_to_matrix(m, orbital_quat);
  glMultMatrixd(m);
  if (orbital_zoom > 2.0) { // gradual pan towards the lander when zoomed in
    sf = 1.0 - exp((2.0-orbital_zoom)/5.0);
    glTranslated(-sf*position.x, -sf*position.y, -sf*position.z);
  }

  if (static_lighting) {
    // Specify light positions here, to fix them in the world coordinate system
    glLightfv(GL_LIGHT2, GL_POSITION, minus_y);
    glLightfv(GL_LIGHT3, GL_POSITION, plus_y);
  }

  // Draw planet
  glColor3f(0.63, 0.33, 0.22);
  glLineWidth(1.0);
  glPushMatrix();
  glRotated(360.0*simulation_time/MARS_DAY, 0.0, 0.0, 1.0); // to make the planet spin
  if (orbital_zoom > 1.0) {
    slices = (int)(16*orbital_zoom); if (slices > 160) slices = 160;
    stacks = (int)(10*orbital_zoom); if (stacks > 100) stacks = 100;
  } else {
    slices = 16; stacks = 10;
  }
  gluQuadricDrawStyle(quadObj, GLU_FILL);
  gluSphere(quadObj, (1.0 - 0.01/orbital_zoom)*MARS_RADIUS, slices, stacks);
  glColor3f(0.31, 0.16, 0.11);
  gluQuadricDrawStyle(quadObj, GLU_LINE);
  gluSphere(quadObj, MARS_RADIUS, slices, stacks);
  glPopMatrix();

  // Draw previous lander positions in cyan that fades with time
  glDisable(GL_LIGHTING);
  glEnable(GL_BLEND);
  glLineWidth(1.0);
  glBegin(GL_LINE_STRIP);
  glColor3f(0.0, 1.0, 1.0);
  glVertex3d(position.x, position.y, position.z);
  j = (track.p+N_TRACK-1)%N_TRACK;
  for (i=0; i<track.n; i++) {
    glColor4f(0.0, 0.75*(N_TRACK-i)/N_TRACK, 0.75*(N_TRACK-i)/N_TRACK, 1.0*(N_TRACK-i)/N_TRACK);
    glVertex3d(track.pos[j].x, track.pos[j].y, track.pos[j].z); 
    j = (j+N_TRACK-1)%N_TRACK;
  }
  glEnd();
  glDisable(GL_BLEND);

  // Draw lander as a cyan dot
  glColor3f(0.0, 1.0, 1.0);
  glPointSize(3.0);
  glBegin(GL_POINTS);
  glVertex3d(position.x, position.y, position.z);
  glEnd();
  glEnable(GL_LIGHTING);

  // Help information
  if (help) display_help_text();

  glutSwapBuffers();
}

void draw_parachute_quad (double d)
  // Draws a single parachute quad, distance d behind the lander
{
  glNormal3d(-1.0, 0.0, 0.0);
  glBegin(GL_QUADS);
  glVertex3d(-d, -LANDER_SIZE, -LANDER_SIZE);
  glVertex3d(-d, -LANDER_SIZE, LANDER_SIZE);
  glVertex3d(-d, LANDER_SIZE, LANDER_SIZE);
  glVertex3d(-d, LANDER_SIZE, -LANDER_SIZE);
  glEnd();
  glBegin(GL_LINES);
  glVertex3d(0.0, 0.0, 0.0);
  glVertex3d(-d, 0.0, 0.0);
  glEnd();
}

void draw_parachute (double d)
  // OpenGL quads and lines to draw a simple parachute, distance d behind the lander
{
  glLineWidth(1.0);
  glColor3f(1.0, 0.75, 0.0);
  glDisable(GL_CULL_FACE);
  draw_parachute_quad(d);
  glPushMatrix();
  glRotated((360.0/M_PI)*atan2(LANDER_SIZE, d), 0.0, 0.0, 1.0);
  draw_parachute_quad(d);
  glPopMatrix();
  glPushMatrix();
  glRotated(-(360.0/M_PI)*atan2(LANDER_SIZE, d), 0.0, 0.0, 1.0);
  draw_parachute_quad(d);
  glPopMatrix();
  glPushMatrix();
  glRotated((360.0/M_PI)*atan2(LANDER_SIZE, d), 0.0, 1.0, 0.0);
  draw_parachute_quad(d);
  glPopMatrix();
  glPushMatrix();
  glRotated(-(360.0/M_PI)*atan2(LANDER_SIZE, d), 0.0, 1.0, 0.0);
  draw_parachute_quad(d);
  glPopMatrix();
  glEnable(GL_CULL_FACE);
}

bool generate_terrain_texture (void)
  // Generates random texture map for surface terrain, with mipmap to avoid aliasing at the horizon
{
  unsigned char *tex_image;
  unsigned long x;
  GLsizei ts;
  bool texture_ok;

  ts = TERRAIN_TEXTURE_SIZE;
  texture_ok = false;
  tex_image = (unsigned char*) calloc(sizeof(unsigned char), TERRAIN_TEXTURE_SIZE*TERRAIN_TEXTURE_SIZE);
  for (x=0; x<TERRAIN_TEXTURE_SIZE*TERRAIN_TEXTURE_SIZE; x++) tex_image[x] = 192 + (unsigned char) (63.0*rand()/RAND_MAX);
  glGenTextures(1, &terrain_texture);
  glBindTexture(GL_TEXTURE_2D, terrain_texture);
  while (!texture_ok && (ts >= 256)) { // try progressively smaller texture maps, give up below 256x256
    glGetError(); // clear error
    if (!gluBuild2DMipmaps(GL_TEXTURE_2D, GL_LUMINANCE, ts, ts, GL_LUMINANCE, GL_UNSIGNED_BYTE, tex_image) && (glGetError() == GL_NO_ERROR)) texture_ok = true;
    else ts /= 2;
  }
  free(tex_image);
  if (texture_ok) {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    return true;
  } else return false;
}

void update_closeup_coords (void)
  // Updates the close-up view's coordinate frame, based on the lander's current position and velocity.
  // This needs to be called every time step, even if the view is not being rendered, since any-angle
  // attitude stabilizers reference closeup_coords.right
{
  vector3d s, tv, t;
  double tmp;

  // Direction from surface to lander (radial) - this must map to the world y-axis
  s = position.norm();

  // Direction of tangential velocity - this must map to the world x-axis
  tv = velocity_from_positions - (velocity_from_positions*s)*s;
  if (tv.abs() < SMALL_NUM) // vertical motion only, use last recorded tangential velocity
    tv = closeup_coords.backwards ? (closeup_coords.right*s)*s - closeup_coords.right : closeup_coords.right - (closeup_coords.right*s)*s; 
  if (tv.abs() > SMALL_NUM) t = tv.norm();

  // Check these two vectors are non-zero and perpendicular (they should be, unless s and closeup_coords.right happen to be parallel)
  if ((tv.abs() <= SMALL_NUM) || (fabs(s*t) > SMALL_NUM)) {
    // Set t to something perpendicular to s
    t.x = -s.y; t.y = s.x; t.z = 0.0;
    if (t.abs() < SMALL_NUM) {t.x = -s.z; t.y = 0.0; t.z = s.x;}
    t = t.norm();
  }

  // Adjust the terrain texture angle if the lander has changed direction. The motion will still be along
  // the x-axis, so we need to rotate the texture to compensate.
  if (closeup_coords.initialized) {
    if (closeup_coords.backwards) {
      tmp = -closeup_coords.right*t;
      if (tmp > 1.0) tmp = 1.0; if (tmp < -1.0) tmp = -1.0;
      if ((-closeup_coords.right^t)*position.norm() < 0.0) terrain_angle += (180.0/M_PI)*acos(tmp);
      else terrain_angle -= (180.0/M_PI)*acos(tmp);
    } else {
      tmp = closeup_coords.right*t;
      if (tmp > 1.0) tmp = 1.0; if (tmp < -1.0) tmp = -1.0;
      if ((closeup_coords.right^t)*position.norm() < 0.0) terrain_angle += (180.0/M_PI)*acos(tmp);
      else terrain_angle -= (180.0/M_PI)*acos(tmp);
    }
    while (terrain_angle < 0.0) terrain_angle += 360.0;
    while (terrain_angle >= 360.0) terrain_angle -= 360.0;
  }

  // Normally we maintain motion to the right, the one exception being when the ground speed passes
  // through zero and changes sign. A sudden 180 degree change of viewpoint would be confusing, so
  // in this instance we allow the lander to fly to the left.
  if (closeup_coords.initialized && (closeup_coords.right*t < 0.0)) {
    closeup_coords.backwards = true;
    closeup_coords.right = -1.0*t;
  } else {
    closeup_coords.backwards = false;
    closeup_coords.right = t;
    closeup_coords.initialized = true;
  }
}

void draw_closeup_window (void)
  // Draws the close-up view of the lander
{
  static double terrain_offset_x = 0.0;
  static double terrain_offset_y = 0.0;
  static double ground_line_offset = 0.0;
  static double last_redraw_time = 0.0;
  static unsigned short rn = 0;
  vector3d s, t, n;
  double lander_drag, chute_drag, glow_factor, aspect_ratio, view_depth, f, tmp, cs, gs;
  double horizon, fog_density, cx, cy, m[16], m2[16], transition_altitude, ground_plane_size;
  unsigned short i, j, rtmp;
  GLfloat fogcolour[4];
  bool dark_side;
  float rand_tri[8];

  glutSetWindow(closeup_window);
  aspect_ratio = (double)view_width/view_height;
  if (do_texture) transition_altitude = TRANSITION_ALTITUDE;
  else transition_altitude = TRANSITION_ALTITUDE_NO_TEXTURE;
  ground_plane_size = 5.0*transition_altitude;

  // Work out an atmospheric haze colour based on prevailing atmospheric density. The power law in the
  // expression below couples with the fog calculation further down, to ensure that the fog doesn't dim
  // the scene on the way down.
  tmp = pow(atmospheric_density(position)/atmospheric_density(vector3d(MARS_RADIUS, 0.0, 0.0)), 0.5);
  if (static_lighting) tmp *= 0.5 * (1.0 + position.norm()*vector3d(0.0, -1.0, 0.0)); // set sky colour
  fogcolour[0] = tmp*0.98; fogcolour[1] = tmp*0.67; fogcolour[2] = tmp*0.52; fogcolour[3] = 0.0;
  glClearColor(tmp*0.98, tmp*0.67, tmp*0.52, 0.0);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (altitude < 0.0) { // just blank the screen if the lander is below the surface
    glutSwapBuffers();
    return;
  }
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  // Set projection matrix and fog values based on altitude.
  // Above the exosphere we see a long way and there is no fog.
  // Between the exosphere and transition_altitude, we see up to the horizon with increasing fog.
  // At transition_altitude we have a totally opaque haze, to disguise the transition from spherical surface to flat surface.
  // Below transition_altitude, we can see as far as the horizon (or transition_altitude with no terrain texture), 
  // with the fog decreasing towards touchdown.
  if (altitude > EXOSPHERE) gluPerspective(CLOSEUP_VIEW_ANGLE, aspect_ratio, 1.0, closeup_offset + 2.0*MARS_RADIUS);
  else {
    horizon = sqrt(position.abs2() - MARS_RADIUS*MARS_RADIUS);
    if (altitude > transition_altitude) {
      f = (altitude-transition_altitude) / (EXOSPHERE-transition_altitude);
      if (f < SMALL_NUM) fog_density = 1000.0; else fog_density = (1.0-f) / (f*horizon);
      view_depth = closeup_offset + horizon;
    } else {
      f = 1.0 - (altitude / transition_altitude);
      if (f < SMALL_NUM) fog_density = 1000.0; else fog_density = (1.0-f) / (f*transition_altitude);
      if (do_texture) {
	fog_density = 0.00005 + 0.5*fog_density;
	view_depth = closeup_offset + horizon;
      } else view_depth = closeup_offset + transition_altitude;
    }
    gluPerspective(CLOSEUP_VIEW_ANGLE, aspect_ratio, 1.0, view_depth);
    glFogf(GL_FOG_DENSITY, fog_density);
    glFogfv(GL_FOG_COLOR, fogcolour);
    if (do_texture) glHint(GL_FOG_HINT, GL_NICEST);
    else glHint(GL_FOG_HINT, GL_FASTEST);
    glEnable(GL_FOG);
  }

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // The world coordinate system in this view is centered on the lander, with y-axis vertical
  // (from planet centre to lander) and x-axis parallel to the lander's tangential velocity.
  // We now need a modelling rotation transformation to map from this system to the planetary
  // coordinate system.

  // Direction from surface to lander (radial) - this must map to the world y-axis
  s = position.norm();

  // Direction of tangential velocity - this must map to the world x-axis
  t = closeup_coords.backwards ? -closeup_coords.right : closeup_coords.right; 

  // Mutual perpendicular to these two vectors - this must map to the world z-axis
  n = t^s;

  // Construct modelling matrix (rotation only) from these three vectors
  m[0] = t.x; m[1] = t.y; m[2] = t.z; m[3] = 0.0;
  m[4] = s.x; m[5] = s.y; m[6] = s.z; m[7] = 0.0;
  m[8] = n.x; m[9] = n.y; m[10] = n.z; m[11] = 0.0;
  m[12] = 0.0; m[13] = 0.0; m[14] = 0.0; m[15] = 1.0;
  invert(m, m2); 

  // Update terrain texture/line offsets
  if (simulation_time != last_redraw_time) {
    terrain_offset_x += cos(terrain_angle*M_PI/180.0) * ground_speed * (simulation_time-last_redraw_time) / (2.0*ground_plane_size);
    terrain_offset_y += sin(terrain_angle*M_PI/180.0) * ground_speed * (simulation_time-last_redraw_time) / (2.0*ground_plane_size);
    while (terrain_offset_x < 0.0) terrain_offset_x += 1.0;
    while (terrain_offset_x >= 1.0) terrain_offset_x -= 1.0;
    while (terrain_offset_y < 0.0) terrain_offset_y += 1.0;
    while (terrain_offset_y >= 1.0) terrain_offset_y -= 1.0;
    if (closeup_coords.backwards) ground_line_offset += ground_speed * (simulation_time-last_redraw_time);
    else ground_line_offset -= ground_speed * (simulation_time-last_redraw_time);
    ground_line_offset -= GROUND_LINE_SPACING*((int)ground_line_offset/(int)(GROUND_LINE_SPACING));
    last_redraw_time = simulation_time;
  }

  // Viewing transformation
  glTranslated(0.0, 0.0, -closeup_offset);
  glRotated(closeup_xr, 1.0, 0.0, 0.0);
  glRotated(closeup_yr + closeup_coords.backwards*180.0, 0.0, 1.0, 0.0);

  if (static_lighting) {
    // Specify light positions here, to fix them in the planetary coordinate system
    glPushMatrix();
    glMultMatrixd(m2); // now in the planetary coordinate system
    glLightfv(GL_LIGHT2, GL_POSITION, minus_y);
    glLightfv(GL_LIGHT3, GL_POSITION, plus_y);
    glLightfv(GL_LIGHT4, GL_POSITION, plus_y);
    glLightfv(GL_LIGHT5, GL_POSITION, plus_z);
    glPopMatrix(); // back to the view's world coordinate system
  }

  // Surface colour
  glColor3f(0.63, 0.33, 0.22);

  if (altitude < transition_altitude) {

    // Draw ground plane below the lander's current position - we need to do this in quarters, with a vertex
    // nearby, to get the fog calculations correct in all OpenGL implementations.
    glBindTexture(GL_TEXTURE_2D, terrain_texture);
    if (do_texture) glEnable(GL_TEXTURE_2D);
    glNormal3d(0.0, 1.0, 0.0);
    glPushMatrix();
    glRotated(terrain_angle, 0.0, 1.0, 0.0);
    glBegin(GL_QUADS);
    glTexCoord2f(1.0 + terrain_offset_x, 1.0 + terrain_offset_y); glVertex3d(ground_plane_size, -altitude, ground_plane_size);      
    glTexCoord2f(1.0 + terrain_offset_x, 0.5 + terrain_offset_y); glVertex3d(ground_plane_size, -altitude, 0.0);
    glTexCoord2f(0.5 + terrain_offset_x, 0.5 + terrain_offset_y); glVertex3d(0.0, -altitude, 0.0);      
    glTexCoord2f(0.5 + terrain_offset_x, 1.0 + terrain_offset_y); glVertex3d(0.0, -altitude, ground_plane_size);
    glTexCoord2f(0.5 + terrain_offset_x, 0.5 + terrain_offset_y); glVertex3d(0.0, -altitude, 0.0);      
    glTexCoord2f(1.0 + terrain_offset_x, 0.5 + terrain_offset_y); glVertex3d(ground_plane_size, -altitude, 0.0);
    glTexCoord2f(1.0 + terrain_offset_x, 0.0 + terrain_offset_y); glVertex3d(ground_plane_size, -altitude, -ground_plane_size);
    glTexCoord2f(0.5 + terrain_offset_x, 0.0 + terrain_offset_y); glVertex3d(0.0, -altitude, -ground_plane_size);
    glTexCoord2f(0.5 + terrain_offset_x, 0.5 + terrain_offset_y); glVertex3d(0.0, -altitude, 0.0);      
    glTexCoord2f(0.5 + terrain_offset_x, 0.0 + terrain_offset_y); glVertex3d(0.0, -altitude, -ground_plane_size);
    glTexCoord2f(0.0 + terrain_offset_x, 0.0 + terrain_offset_y); glVertex3d(-ground_plane_size, -altitude, -ground_plane_size);
    glTexCoord2f(0.0 + terrain_offset_x, 0.5 + terrain_offset_y); glVertex3d(-ground_plane_size, -altitude, 0.0);
    glTexCoord2f(0.5 + terrain_offset_x, 1.0 + terrain_offset_y); glVertex3d(0.0, -altitude, ground_plane_size);
    glTexCoord2f(0.5 + terrain_offset_x, 0.5 + terrain_offset_y); glVertex3d(0.0, -altitude, 0.0);      
    glTexCoord2f(0.0 + terrain_offset_x, 0.5 + terrain_offset_y); glVertex3d(-ground_plane_size, -altitude, 0.0);
    glTexCoord2f(0.0 + terrain_offset_x, 1.0 + terrain_offset_y); glVertex3d(-ground_plane_size, -altitude, ground_plane_size);
    glEnd();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_DEPTH_TEST);

    if (!do_texture) { // draw lines on the ground plane at constant x (to show ground speed)
      glEnable(GL_BLEND);
      glLineWidth(2.0);
      glBegin(GL_LINES);
      if (closeup_coords.backwards) tmp = -ground_line_offset - transition_altitude;
      else tmp = ground_line_offset + transition_altitude;
      while ((closeup_coords.backwards ? -tmp : tmp) > -transition_altitude) {
	// Fade the lines out towards the horizon, to avoid aliasing artefacts. The fade is a function of distance from the
	// centre (tmp) and altitude: the lower the lander gets, the more pronounced the fade.
	// We need to do draw each line in two parts, with a vertex nearby, to get the fog calculations correct in all OpenGL implementations.
	// To make the lines fade more strongly when landed, decrease the second number.
	// To make the lines less apparent at high altitude, decrease the first number. 
	f = exp( -fabs( pow((transition_altitude-altitude) / transition_altitude, 10.0) * tmp / (10.0*GROUND_LINE_SPACING)) );
	glColor4f(0.32, 0.17, 0.11, f);
	glVertex3d(tmp, -altitude, -transition_altitude);
	glVertex3d(tmp, -altitude, 0.0);
	glVertex3d(tmp, -altitude, 0.0);
	glVertex3d(tmp, -altitude, transition_altitude);
	if (closeup_coords.backwards) tmp += GROUND_LINE_SPACING;
	else tmp -= GROUND_LINE_SPACING;
      }
      glEnd();
      glDisable(GL_BLEND);
    }

    if (!crashed) { // draw a circular shadow below the lander
      glColor3f(0.32, 0.17, 0.11);
      glBegin(GL_TRIANGLES);
      for (i=0; i<360; i+=10) {
	glVertex3d(0.0, -altitude, 0.0);
	glVertex3d(LANDER_SIZE*cos(M_PI*(i+10)/180.0), -altitude, LANDER_SIZE*sin(M_PI*(i+10)/180.0));
	glVertex3d(LANDER_SIZE*cos(M_PI*i/180.0), -altitude, LANDER_SIZE*sin(M_PI*i/180.0));
      }
      glEnd();
    } else {
      rtmp = 0;
      glColor3f(1.0, 1.0, 1.0);
      glBegin(GL_TRIANGLES); // draw some shards of metal
      for (i=0; i<60; i++) {
	for (j=0; j<8; j++) { rand_tri[j] = randtab[rtmp]; rtmp = (rtmp+1)%N_RAND; }
	cx = 40.0 * (rand_tri[0] - 0.5);
	cy = 40.0 * (rand_tri[1] - 0.5);
	glNormal3d(0.0, 1.0, 0.0);
	glVertex3d(cx + 2.0*LANDER_SIZE*rand_tri[2], -altitude, cy + 2.0*LANDER_SIZE*rand_tri[3]);
	glVertex3d(cx + 2.0*LANDER_SIZE*rand_tri[4], -altitude, cy + 2.0*LANDER_SIZE*rand_tri[5]);
	glVertex3d(cx + 2.0*LANDER_SIZE*rand_tri[6], -altitude, cy + 2.0*LANDER_SIZE*rand_tri[7]);
      }
      glEnd();
      if (parachute_status != LOST) {
	glColor3f(1.0, 1.0, 0.0);
	glBegin(GL_TRIANGLES);  // draw some shreds of yellow canvas
	for (i=0; i<30; i++) {
	  for (j=0; j<8; j++) { rand_tri[j] = randtab[rtmp]; rtmp = (rtmp+1)%N_RAND; }
	  cx = 40.0 * (rand_tri[0] - 0.5);
	  cy = 40.0 * (rand_tri[1] - 0.5);
	  glNormal3d(0.0, 1.0, 0.0);
	  glVertex3d(cx + 2.0*LANDER_SIZE*rand_tri[2], -altitude, cy + 2.0*LANDER_SIZE*rand_tri[3]);
	  glVertex3d(cx + 2.0*LANDER_SIZE*rand_tri[4], -altitude, cy + 2.0*LANDER_SIZE*rand_tri[5]);
	  glVertex3d(cx + 2.0*LANDER_SIZE*rand_tri[6], -altitude, cy + 2.0*LANDER_SIZE*rand_tri[7]);
	}
	glEnd();
      }
    }
    glEnable(GL_DEPTH_TEST);
  
  } else {

    // Draw spherical planet - can disable depth test (for speed)
    glDisable(GL_DEPTH_TEST);
    glPushMatrix();

    if (altitude > EXOSPHERE) {

      // Draw the planet reduced size at a reduced displacement, to avoid numerical OpenGL problems with huge viewing distances.
      glTranslated(0.0, -MARS_RADIUS, 0.0);
      glMultMatrixd(m2); // now in the planetary coordinate system
      glRotated(360.0*simulation_time/MARS_DAY, 0.0, 0.0, 1.0); // to make the planet spin
      glutMottledSphere(MARS_RADIUS * (MARS_RADIUS / (altitude + MARS_RADIUS)), 160, 100);

    } else {

      // Draw the planet actual size at the correct displacement
      glTranslated(0.0, -(MARS_RADIUS + altitude), 0.0);
      glMultMatrixd(m2); // now in the planetary coordinate system
      glRotated(360.0*simulation_time/MARS_DAY, 0.0, 0.0, 1.0); // to make the planet spin
      glutMottledSphere(MARS_RADIUS, 160, 100);

    }

    glPopMatrix(); // back to the view's world coordinate system
    glEnable(GL_DEPTH_TEST);

  }

  glDisable(GL_FOG); // fog only applies to the ground
  dark_side = (static_lighting && (position.y > 0.0) && (sqrt(position.x*position.x + position.z*position.z) < MARS_RADIUS));
  if (dark_side) { // in the shadow of the planet, we need some diffuse lighting to highlight the lander
    glDisable(GL_LIGHT2); glDisable(GL_LIGHT3); 
    glEnable(GL_LIGHT4); glEnable(GL_LIGHT5);
  }

  // Work out drag on lander - if it's high, we will surround the lander with an incandescent glow. Also
  // work out drag on parachute: if it's zero, we will not draw the parachute fully open behind the lander.
  // Assume high Reynolds number, quadratic drag = -0.5 * rho * v^2 * A * C_d
  lander_drag = 0.5*DRAG_COEF_LANDER*atmospheric_density(position)*M_PI*LANDER_SIZE*LANDER_SIZE*velocity_from_positions.abs2();
  chute_drag = 0.5*DRAG_COEF_CHUTE*atmospheric_density(position)*5.0*2.0*LANDER_SIZE*2.0*LANDER_SIZE*velocity_from_positions.abs2();

  // Draw the lander's parachute - behind the lander in the direction of travel
  if ( (parachute_status == DEPLOYED) && !crashed ) {
    if (velocity_from_positions.abs() < SMALL_NUM) {
      // Lander is apparently stationary - so draw the parachute above and near to the lander
      gs = 0.0; cs = -1.0; tmp = 2.0;
    } else {
      gs = ground_speed; cs = climb_speed;
      if (chute_drag) tmp = 5.0; // parachute fully open
      else tmp = 2.0; // parachute not fully open
    }
    glPushMatrix();
    glRotated((180.0/M_PI)*atan2(cs, gs), 0.0, 0.0, 1.0);
    draw_parachute(tmp);
    glPopMatrix();
  }

  // Display help arrow to show surface direction
  if (help) display_help_arrows();
  else display_help_prompt();

  // Switch to the planetary coordinate system
  glPushMatrix();
  glMultMatrixd(m2);

  // Lander orientation relative to planetary coordinate system - xyz Euler angles
  xyz_euler_to_matrix(orientation, m);
  glMultMatrixd(m);

  // Put lander's centre of gravity at the origin
  glTranslated(0.0, 0.0, -LANDER_SIZE/2);

  // Draw lander
  if (!crashed) {
    glColor3f(1.0, 1.0, 1.0);
    glutCone(LANDER_SIZE, LANDER_SIZE, 50, 50, true);
  }

  if (dark_side) { // back to standard lighting model
    glEnable(GL_LIGHT2); glEnable(GL_LIGHT3);
    glDisable(GL_LIGHT4); glDisable(GL_LIGHT5);
  }

  // Draw engine exhaust flare
  if (thrust_wrt_world().abs() > 0.0) {
    glColor3f(1.0, 0.5, 0.0);
    glRotated(180.0, 1.0, 0.0, 0.0);
    glDisable(GL_LIGHTING);
    glutCone(LANDER_SIZE/2, 2*LANDER_SIZE*thrust_wrt_world().abs()/MAX_THRUST, 50, 50, false);
    glEnable(GL_LIGHTING);
  }

  glPopMatrix(); // back to the world coordinate system

  // Draw incandescent glow surrounding lander
  if (lander_drag*velocity_from_positions.abs() > HEAT_FLUX_GLOW_THRESHOLD) {
    // Calculate an heuristic "glow factor", in the range 0 to 1, for graphics effects
    glow_factor = (lander_drag*velocity_from_positions.abs()-HEAT_FLUX_GLOW_THRESHOLD) / (4.0*HEAT_FLUX_GLOW_THRESHOLD); 
    if (glow_factor > 1.0) glow_factor = 1.0;
    glow_factor *= 0.7 + 0.3*randtab[rn]; rn = (rn+1)%N_RAND; // a little random variation for added realism
    glRotated((180.0/M_PI)*atan2(climb_speed, ground_speed), 0.0, 0.0, 1.0);
    glRotated(-90.0, 0.0, 1.0, 0.0);
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glColor4f(1.0, glow_factor, 0.0, 0.8*glow_factor);
    glutCone(1.25*LANDER_SIZE, (2.0 + 10.0*glow_factor)*LANDER_SIZE, 50, 50+(int)(250*glow_factor), false);
    glutOpenHemisphere(1.25*LANDER_SIZE, 50, 50);
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
  }

  glutSwapBuffers(); 
}

void draw_main_window (void)
  // Draw grey lines to partition the display into three sub-windows
{
  glutSetWindow(main_window);
  glClear(GL_COLOR_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glColor3f(0.5, 0.5, 0.5);
  glBegin(GL_LINES);
  glVertex2i(view_width + 2*GAP, INSTRUMENT_HEIGHT + 2*GAP);
  glVertex2i(view_width + 2*GAP, win_height);
  glVertex2i(0, INSTRUMENT_HEIGHT + 2*GAP);
  glVertex2i(win_width, INSTRUMENT_HEIGHT + 2*GAP);
  glEnd();
  glutSwapBuffers();
}

void refresh_all_subwindows (void)
  // Marks all subwindows as needing a redraw every n times called, where n depends on the simulation speed
{
  static unsigned short n = 0;

  if (simulation_speed > 5) {
    switch (simulation_speed) {
    case 6:
      n += 500; // update graphics every 2 iterations
      break;
    case 7:
      n += 100; // update graphics every 10 iterations
      break;
    case 8:
      n += 20; // update graphics every 50 iterations
      break;
    case 9:
      n += 10;  // update graphics every 100 iterations
      break;
    case 10:
      n += 1;  // update graphics every 1000 iterations
      break;
    }
    if (n>=1000) n=0;
    if (!paused && !landed && n) return;
  }

  glutPostWindowRedisplay(closeup_window);
  glutPostWindowRedisplay(orbital_window);
  glutPostWindowRedisplay(instrument_window);
}

bool safe_to_deploy_parachute (void)
  // Checks whether the parachute is safe to deploy at the current position and velocity
{
  double drag;

  // Assume high Reynolds number, quadratic drag = -0.5 * rho * v^2 * A * C_d
  drag = 0.5*DRAG_COEF_CHUTE*atmospheric_density(position)*5.0*2.0*LANDER_SIZE*2.0*LANDER_SIZE*velocity_from_positions.abs2();
  // Do not use the global variable "altitude" here, in case this function is called from within the
  // numerical_dynamics function, before altitude is updated in the update_visualization function
  if ((drag > MAX_PARACHUTE_DRAG) || ((velocity_from_positions.abs() > MAX_PARACHUTE_SPEED) && ((position.abs() - MARS_RADIUS) < EXOSPHERE))) return false;
  else return true;
}

void update_visualization (void)
  // The visualization part of the idle function. Re-estimates altitude, velocity, climb speed and ground
  // speed from current and previous positions. Updates throttle and fuel levels, then redraws all subwindows.
{
  static vector3d last_track_position;
  vector3d av_p, d;
  double a, b, c, mu;

  simulation_time += delta_t;
  altitude = position.abs() - MARS_RADIUS;

  // Use average of current and previous positions when calculating climb and ground speeds
  av_p = (position + last_position).norm();
  if (delta_t != 0.0) velocity_from_positions = (position - last_position)/delta_t;
  else velocity_from_positions = vector3d(0.0, 0.0, 0.0);
  climb_speed = velocity_from_positions*av_p;
  ground_speed = (velocity_from_positions - climb_speed*av_p).abs();

  // Check to see whether the lander has landed
  if (altitude < LANDER_SIZE/2.0) {
    glutIdleFunc(NULL);
    // Estimate position and time of impact
    d = position - last_position;
    a = d.abs2();
    b = 2.0*last_position*d;
    c = last_position.abs2() - (MARS_RADIUS + LANDER_SIZE/2.0) * (MARS_RADIUS + LANDER_SIZE/2.0);
    mu = (-b - sqrt(b*b-4.0*a*c))/(2.0*a);
    position = last_position + mu*d;
    simulation_time -= (1.0-mu)*delta_t; 
    altitude = LANDER_SIZE/2.0;
    landed = true;
    if ((fabs(climb_speed) > MAX_IMPACT_DESCENT_RATE) || (fabs(ground_speed) > MAX_IMPACT_GROUND_SPEED)) crashed = true;
    velocity_from_positions = vector3d(0.0, 0.0, 0.0);
  }

  // Update throttle and fuel (throttle might have been adjusted by the autopilot)
  if (throttle < 0.0) throttle = 0.0;
  if (throttle > 1.0) throttle = 1.0;
  fuel -= delta_t * (FUEL_RATE_AT_MAX_THRUST*throttle) / FUEL_CAPACITY;
  if (fuel <= 0.0) fuel = 0.0;
  if (landed || (fuel == 0.0)) throttle = 0.0;
  throttle_control = (short)(throttle*THROTTLE_GRANULARITY + 0.5);

  // Check to see whether the parachute has vaporized or the tethers have snapped
  if (parachute_status == DEPLOYED) {
    if (!safe_to_deploy_parachute() || parachute_lost) {
      parachute_lost = true; // to guard against the autopilot reinstating the parachute!
      parachute_status = LOST;
    }
  }

  // Update record of lander's previous positions, but only if the position or the velocity has 
  // changed significantly since the last update
  if ( !track.n || (position-last_track_position).norm() * velocity_from_positions.norm() < TRACK_ANGLE_DELTA
      || (position-last_track_position).abs() > TRACK_DISTANCE_DELTA ) {
    track.pos[track.p] = position;
    track.n++; if (track.n > N_TRACK) track.n = N_TRACK;
    track.p++; if (track.p == N_TRACK) track.p = 0;
    last_track_position = position;
  }

  // Redraw everything
  refresh_all_subwindows();
}

void attitude_stabilization (void)
  // Three-axis stabilization to ensure the lander's base is always pointing downwards 
{
  vector3d up, left, out;
  double m[16];

  up = position.norm(); // this is the direction we want the lander's nose to point in

  // !!!!!!!!!!!!! HINT TO STUDENTS ATTEMPTING THE EXTENSION EXERCISES !!!!!!!!!!!!!!
  // For any-angle attitude control, we just need to set "up" to something different,
  // and leave the remainder of this function unchanged. For example, suppose we want
  // the attitude to be stabilized at stabilized_attitude_angle to the vertical in the
  // close-up view. So we need to rotate "up" by stabilized_attitude_angle degrees around
  // an axis perpendicular to the plane of the close-up view. This axis is given by the
  // vector product of "up"and "closeup_coords.right". To calculate the result of the
  // rotation, search the internet for information on the axis-angle rotation formula.

  //In plane attitude control
  vector3d k = (up ^ closeup_coords.right).norm();
  double saa_rad = stabilized_attitude_angle * M_PI / 180;
  up = up*cos(saa_rad) + (k^up)*sin(saa_rad) + k*(k*up)*(1-cos(saa_rad));
  closeup_coords_right = closeup_coords.right;

  //Complete attitude control, rotate up about position.norm()
  vector3d k2 = position.norm();//(k ^ position).norm();
  double oopa_rad = out_of_plane_angle * M_PI / 180;
  up = up*cos(oopa_rad) + (k2 ^ up)*sin(oopa_rad) + k2*(k2*up)*(1 - cos(oopa_rad));

  // Set left to something perpendicular to up
  left.x = -up.y; left.y = up.x; left.z = 0.0;
  if (left.abs() < SMALL_NUM) {left.x = -up.z; left.y = 0.0; left.z = up.x;}
  left = left.norm();  
  out = left^up;
  // Construct modelling matrix (rotation only) from these three vectors
  m[0] = out.x; m[1] = out.y; m[2] = out.z; m[3] = 0.0;
  m[4] = left.x; m[5] = left.y; m[6] = left.z; m[7] = 0.0;
  m[8] = up.x; m[9] = up.y; m[10] = up.z; m[11] = 0.0;
  m[12] = 0.0; m[13] = 0.0; m[14] = 0.0; m[15] = 1.0;
  // Decomponse into xyz Euler angles
  orientation = matrix_to_xyz_euler(m);
}

vector3d thrust_wrt_world (void)
  // Works out thrust vector in the world reference frame, given the lander's orientation
{
  double m[16], k, delayed_throttle, lag = ENGINE_LAG;
  vector3d a, b;
  static double lagged_throttle = 0.0;
  static double last_time_lag_updated = -1.0;

  if (simulation_time < last_time_lag_updated) lagged_throttle = 0.0; // simulation restarted
  if (throttle < 0.0) throttle = 0.0;
  if (throttle > 1.0) throttle = 1.0;
  if (landed || (fuel == 0.0)) throttle = 0.0;

  if (simulation_time != last_time_lag_updated) {

    // Delayed throttle value from the throttle history buffer
    if (throttle_buffer_length > 0) {
      delayed_throttle = throttle_buffer[throttle_buffer_pointer];
      throttle_buffer[throttle_buffer_pointer] = throttle;
      throttle_buffer_pointer = (throttle_buffer_pointer + 1) % throttle_buffer_length;
    } else delayed_throttle = throttle;

    // Lag, with time constant ENGINE_LAG
    if (lag <= 0.0) k = 0.0;
    else k = pow(exp(-1.0), delta_t/lag);
    lagged_throttle = k*lagged_throttle + (1.0-k)*delayed_throttle;

    last_time_lag_updated = simulation_time;
  }

  if (stabilized_attitude && (stabilized_attitude_angle == 0)) { // specific solution, avoids rounding errors in the more general calculation below
    b = lagged_throttle*MAX_THRUST*position.norm();
  } else {
    a.x = 0.0; a.y = 0.0; a.z = lagged_throttle*MAX_THRUST;
    xyz_euler_to_matrix(orientation, m);
    b.x = m[0]*a.x + m[4]*a.y + m[8]*a.z;
    b.y = m[1]*a.x + m[5]*a.y + m[9]*a.z;
    b.z = m[2]*a.x + m[6]*a.y + m[10]*a.z;
  }
  return b;
}

void update_lander_state (void)
  // The GLUT idle function, called every time round the event loop
{
  unsigned long delay;

  // User-controlled delay
  if ((simulation_speed > 0) && (simulation_speed < 5)) {
    delay = (5-simulation_speed)*MAX_DELAY/4;
#ifdef WIN32
    Sleep(delay/1000); // milliseconds
#else
    usleep( (useconds_t)delay ); // microseconds
#endif
  }

  // This needs to be called every time step, even if the close-up view is not being rendered,
  // since any-angle attitude stabilizers reference closeup_coords.right
  update_closeup_coords();
  
  // Update historical record
  last_position = position;

  // Mechanical dynamics
  numerical_dynamics();

  // Refresh the visualization
  update_visualization();
}

void reset_simulation (void)
  // Resets the simulation to the initial state
{
  vector3d p, tv;
  unsigned long i;

  // Reset these three lander parameters here, so they can be overwritten in initialize_simulation() if so desired
  stabilized_attitude_angle = 0;
  out_of_plane_angle = 0;
  throttle = 0.0;
  fuel = 1.0;

  // Restore initial lander state
  initialize_simulation();

  // Check whether the lander is underground - if so, make sure it doesn't move anywhere
  landed = false;
  crashed = false;
  altitude = position.abs() - MARS_RADIUS;
  if (altitude < LANDER_SIZE/2.0) {
    glutIdleFunc(NULL);
    landed = true;
    velocity = vector3d(0.0, 0.0, 0.0);
  }

  // Visualisation routine's record of various speeds and velocities
  velocity_from_positions = velocity;
  last_position = position - delta_t*velocity_from_positions;
  p = position.norm();
  climb_speed = velocity_from_positions*p;
  tv = velocity_from_positions - climb_speed*p;
  ground_speed = tv.abs();

  // Miscellaneous state variables
  throttle_control = (short)(throttle*THROTTLE_GRANULARITY + 0.5);
  simulation_time = 0.0;  
  track.n = 0;
  parachute_lost = false;
  closeup_coords.initialized = false;
  closeup_coords.backwards = false;
  closeup_coords.right = vector3d(1.0, 0.0, 0.0);
  update_closeup_coords();

  // Initialize the throttle history buffer
  if (delta_t > 0.0) throttle_buffer_length = (unsigned long) (ENGINE_DELAY/delta_t + 0.5);
  else throttle_buffer_length = 0;
  if (throttle_buffer_length > 0) {
    if (throttle_buffer != NULL) delete[] throttle_buffer;
    throttle_buffer = new double[throttle_buffer_length];
    for (i=0; i<throttle_buffer_length; i++) throttle_buffer[i] = throttle;
    throttle_buffer_pointer = 0;
  }

  // Reset GLUT state
  if (paused || landed) refresh_all_subwindows();
  else glutIdleFunc(update_lander_state);
}

void set_orbital_projection_matrix (void)
  // Called from reshape and zoom functions
{
  double aspect_ratio;

  aspect_ratio = (double)view_width/(double)view_height;
  glutSetWindow(orbital_window);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-2.0*MARS_RADIUS*aspect_ratio/orbital_zoom, 2.0*MARS_RADIUS*aspect_ratio/orbital_zoom, 
	  -2.0*MARS_RADIUS/orbital_zoom, 2.0*MARS_RADIUS/orbital_zoom, -100.0*MARS_RADIUS, 100.0*MARS_RADIUS);
}

void reshape_main_window (int width, int height)
  // Called when the main window is created or resized
{
  // Resize the main window
  glutSetWindow(main_window);
  win_width = glutGet(GLUT_WINDOW_WIDTH);
  win_height = glutGet(GLUT_WINDOW_HEIGHT);

  // Work out subwindow dimensions and set projection matrix for main window
  view_width = (win_width - 4*GAP)/2;
  if (view_width < 1) view_width = 1; // GLUT warns about non-positive window dimensions
  view_height = (win_height - INSTRUMENT_HEIGHT - 4*GAP);
  if (view_height < 1) view_height = 1; // GLUT warns about non-positive window dimensions
  glViewport(0, 0, win_width, win_height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, win_width, 0, win_height, -1.0, 1.0);

  // Resize and initialize the close-up view window
  glutSetWindow(closeup_window);
  glutPositionWindow(GAP, GAP);
  glutReshapeWindow(view_width, view_height);
  glViewport(0, 0, view_width, view_height);
  glDrawBuffer(GL_BACK);
  glutPostRedisplay();
  
  // Resize and initialize the orbital view window
  glutSetWindow(orbital_window);
  glutPositionWindow(view_width + 3*GAP, GAP);
  glutReshapeWindow(view_width, view_height);
  glViewport(0, 0, view_width, view_height);
  set_orbital_projection_matrix();
  glDrawBuffer(GL_BACK);
  glutPostRedisplay();

  // Resize and initialize the instrument window
  glutSetWindow(instrument_window);
  glutPositionWindow(GAP, view_height + 3*GAP);
  glutReshapeWindow(2*(view_width+GAP), INSTRUMENT_HEIGHT);
  glViewport(0, 0, 2*(view_width+GAP), INSTRUMENT_HEIGHT);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, 2*(view_width+GAP), 0, INSTRUMENT_HEIGHT, -1.0, 1.0); 
  glDrawBuffer(GL_BACK); 
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glutPostRedisplay();
}

void orbital_mouse_button (int button, int state, int x, int y)
  // Callback for mouse button presses in the orbital view window
{
  if ((button == GLUT_WHEEL_UP) || (((button == GLUT_MIDDLE_BUTTON) || glutGetModifiers()) && (state == GLUT_DOWN))) {
    if (orbital_zoom < 100.0) { // don't let them get too close!
      orbital_zoom /= 0.9;
      if ((button == GLUT_MIDDLE_BUTTON) || glutGetModifiers()) orbital_zoom /= 0.9; // to match wheel events
    }
    save_orbital_zoom = -1.0;
    set_orbital_projection_matrix();
    if (paused || landed) refresh_all_subwindows();
  }
  else if ((button == GLUT_WHEEL_DOWN) || ((button == GLUT_RIGHT_BUTTON) && (state == GLUT_DOWN))) {
    if (orbital_zoom > 0.001) {
      orbital_zoom *= 0.9;
      if (button == GLUT_RIGHT_BUTTON) orbital_zoom *= 0.9; // to match wheel events
    }
    save_orbital_zoom = -1.0;
    set_orbital_projection_matrix();
    if (paused || landed) refresh_all_subwindows();
  }
  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_UP) {
      last_click_x = -1;
      last_click_y = -1;
    }
    if (state == GLUT_DOWN) {
      last_click_x = x;
      last_click_y = y;
    }
  }
}

void orbital_mouse_motion (int x, int y)
  // Callback for mouse drags in the orbital view window
{
  quat_t spin_quat;
  
  if (last_click_x < 0) return; // not the left mouse button

  spin_quat = track_quats((2.0*last_click_x - view_width) / view_width,
			  (view_height - 2.0*last_click_y) / view_height,
			  (2.0*x - view_width) / view_width,
			  (view_height - 2.0*y) / view_height);
  orbital_quat = add_quats(spin_quat, orbital_quat);
  last_click_x = x;
  last_click_y = y;
  if (paused || landed) refresh_all_subwindows();
}

void closeup_mouse_button (int button, int state, int x, int y)
  // Callback for mouse button presses in the close-up view window
{
  if ((button == GLUT_WHEEL_UP) || (((button == GLUT_MIDDLE_BUTTON) || glutGetModifiers()) && (state == GLUT_DOWN))) {
    if (closeup_offset > 2.0*LANDER_SIZE) {
      closeup_offset *= 0.9;
      if ((button == GLUT_MIDDLE_BUTTON) || glutGetModifiers()) closeup_offset *= 0.9; // to match wheel events
    }
    if (paused || landed) refresh_all_subwindows();
  }
  else if ((button == GLUT_WHEEL_DOWN) || ((button == GLUT_RIGHT_BUTTON) && (state == GLUT_DOWN))) {
    if (closeup_offset < 200.0*LANDER_SIZE) {
      closeup_offset /= 0.9;
      if (button == GLUT_RIGHT_BUTTON) closeup_offset /= 0.9; // to match wheel events
    }
    if (paused || landed) refresh_all_subwindows();
  }
  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_UP) {
      last_click_y = -1;
      last_click_x = -1;
    }
    if (state == GLUT_DOWN) {
      last_click_y = y;
      last_click_x = x;
    }
  }
}

void closeup_mouse_motion (int x, int y)
  // Callback for mouse drags in the close-up view window
{
  if (last_click_x < 0) return; // not the left mouse button

  closeup_xr += y - last_click_y;
  closeup_yr += x - last_click_x;
  if (closeup_xr < 0.0) closeup_xr = 0.0;
  if (closeup_xr > 90.0) closeup_xr = 90.0;
  last_click_y = y;
  last_click_x = x;
  if (paused || landed) refresh_all_subwindows();
}

void glut_special (int key, int x, int y)
  // Callback for special key presses in all windows
{
  switch(key) {
  case GLUT_KEY_UP: // throttle up
    if (!autopilot_enabled && !landed && (fuel>0.0)) {
      throttle_control++;
      if (throttle_control>THROTTLE_GRANULARITY) throttle_control = THROTTLE_GRANULARITY;
      throttle = (double)throttle_control/THROTTLE_GRANULARITY;
    }
    break;
  case GLUT_KEY_DOWN: // throttle down
    if (!autopilot_enabled && !landed) {
      throttle_control--;
      if (throttle_control<0) throttle_control = 0;
      throttle = (double)throttle_control/THROTTLE_GRANULARITY;
    }
    break;
  case GLUT_KEY_RIGHT: // faster simulation
    simulation_speed++;
    if (simulation_speed>10) simulation_speed = 10;
    if (paused) {
      if (!landed) glutIdleFunc(update_lander_state);
      paused = false;
    }
    break;
  case GLUT_KEY_LEFT: // slower simulation
    simulation_speed--;
    if (simulation_speed<0) simulation_speed = 0;
    if (!simulation_speed) {
      glutIdleFunc(NULL);
      paused = true;
    }
    break;
  }
  if (paused || landed) refresh_all_subwindows();
}

void glut_key (unsigned char k, int x, int y)
  // Callback for key presses in all windows
{
  switch(k) {

  case 'z':
	  if (!autopilot_enabled && !landed && stabilized_attitude){
		  stabilized_attitude_angle++;
	  }
	  break;

  case 'x':
	  if (!autopilot_enabled && !landed && stabilized_attitude){
		  stabilized_attitude_angle--;
	  }
	  break;

  case 'c':
	  if (!autopilot_enabled && !landed && stabilized_attitude){
		  out_of_plane_angle++;
	  }
	  break;

  case 'v':
	  if (!autopilot_enabled && !landed && stabilized_attitude){
		  out_of_plane_angle--;
	  }
	  break;
    
  case 27: case 'q': case 'Q':
    // Escape or q or Q  - exit
    exit(0);
    break;

  case '0':
    // switch to scenario 0
    scenario = 0;
    reset_simulation();
    break;

  case '1':
    // switch to scenario 1
    scenario = 1;
    reset_simulation();
    break;

  case '2':
    // switch to scenario 2
    scenario = 2;
    reset_simulation();
    break;

  case '3':
    // switch to scenario 3
    scenario = 3;
    reset_simulation();
    break;

  case '4':
    // switch to scenario 4
    scenario = 4;
    reset_simulation();
    break;

  case '5':
    // switch to scenario 5
    scenario = 5;
    reset_simulation();
    break;

  case '6':
    // switch to scenario 6
    scenario = 6;
    reset_simulation();
    break;

  case '7':
    // switch to scenario 7
    scenario = 7;
    reset_simulation();
    break;

  case '8':
    // switch to scenario 8
    scenario = 8;
    reset_simulation();
    break;

  case '9':
    // switch to scenario 9
    scenario = 9;
    reset_simulation();
    break;

  case 'a': case 'A':
    // a or A - autopilot
    if (!landed) autopilot_enabled = !autopilot_enabled;
    if (paused) refresh_all_subwindows();
    break;

  case 'h': case 'H':
    // h or H - help
    if (help) {
      help = false;
      if (save_orbital_zoom > 0.0) orbital_zoom = save_orbital_zoom;
    } else {
      help = true;
      save_orbital_zoom = orbital_zoom;
      orbital_zoom = 0.4;
    }
    set_orbital_projection_matrix();
    if (paused || landed) refresh_all_subwindows();
    break;

  case 'l': case 'L':
    // l or L - toggle lighting model
    static_lighting = !static_lighting;
    glutSetWindow(orbital_window); enable_lights();
    glutSetWindow(closeup_window); enable_lights();
    if (paused || landed) refresh_all_subwindows();
    break;

  case 't': case 'T':
    // t or T - terrain texture
    do_texture = !do_texture;
    if (!texture_available) do_texture = false;
    if (paused || landed) refresh_all_subwindows();
    break;

  case 'p': case 'P':
    // p or P - deploy parachute
    if (!autopilot_enabled && !landed && (parachute_status == NOT_DEPLOYED)) parachute_status = DEPLOYED;
    if (paused) refresh_all_subwindows();
    break;

  case 's': case 'S':
    // s or S - attitude stabilizer
    if (!autopilot_enabled && !landed) stabilized_attitude = !stabilized_attitude;
    if (paused) refresh_all_subwindows();
    break;

  case 32:
    // space bar
    simulation_speed = 0;
    glutIdleFunc(NULL);
    if (paused && !landed) update_lander_state();
    else refresh_all_subwindows();
    paused = true;
    break;

  }
}

int main (int argc, char* argv[])
  // Initializes GLUT windows and lander state, then enters GLUT main loop
{
  int i;

  // Main GLUT window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(PREFERRED_WIDTH, PREFERRED_HEIGHT);
  view_width = (PREFERRED_WIDTH - 4*GAP)/2;
  view_height = (PREFERRED_HEIGHT - INSTRUMENT_HEIGHT - 4*GAP);
  main_window = glutCreateWindow("Mars Lander (Gabor Csanyi and Andrew Gee, October 2014)");
  glDrawBuffer(GL_BACK);
  glLineWidth(2.0);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glutDisplayFunc(draw_main_window);
  glutReshapeFunc(reshape_main_window);  
  glutIdleFunc(update_lander_state);
  glutKeyboardFunc(glut_key);
  glutSpecialFunc(glut_special);

  // The close-up view subwindow
  closeup_window = glutCreateSubWindow(main_window, GAP, GAP, view_width, view_height);
  glDrawBuffer(GL_BACK);
  setup_lights();
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glEnable(GL_CULL_FACE); // we only need back faces for the parachute
  glDisable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_NORMALIZE);
  glDepthFunc(GL_LEQUAL);
  glShadeModel(GL_SMOOTH);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE); // we need two-sided lighting for the parachute
  glEnable(GL_COLOR_MATERIAL);
  glFogi(GL_FOG_MODE, GL_EXP);
  glutDisplayFunc(draw_closeup_window);
  glutMouseFunc(closeup_mouse_button);
  glutMotionFunc(closeup_mouse_motion);
  glutKeyboardFunc(glut_key);
  glutSpecialFunc(glut_special);
  texture_available = generate_terrain_texture();
  if (!texture_available) do_texture = false;
  closeup_offset = 50.0;
  closeup_xr = 10.0;
  closeup_yr = 0.0;
  terrain_angle = 0.0;

  // The orbital view subwindow
  orbital_window = glutCreateSubWindow(main_window, view_width + 3*GAP, GAP, view_width, view_height);
  glDrawBuffer(GL_BACK);
  setup_lights();
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glEnable(GL_CULL_FACE); // since the only polygons in this view define a solid sphere
  glDisable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_NORMALIZE);
  glDepthFunc(GL_LEQUAL);
  glShadeModel(GL_SMOOTH);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glutDisplayFunc(draw_orbital_window);
  glutMouseFunc(orbital_mouse_button);
  glutMotionFunc(orbital_mouse_motion);
  glutKeyboardFunc(glut_key);
  glutSpecialFunc(glut_special);
  quadObj = gluNewQuadric();
  orbital_quat.v.x = 0.53; orbital_quat.v.y = -0.21;
  orbital_quat.v.z = 0.047; orbital_quat.s = 0.82;
  normalize_quat(orbital_quat);
  save_orbital_zoom = 1.0;
  orbital_zoom = 1.0;

  // The instrument subwindow
  instrument_window = glutCreateSubWindow(main_window, GAP, view_height + 3*GAP, 2*(view_width+GAP), INSTRUMENT_HEIGHT);
  glutDisplayFunc(draw_instrument_window);
  glutKeyboardFunc(glut_key);
  glutSpecialFunc(glut_special);

  // Generate the random number table
  srand(0);
  for (i=0; i<N_RAND; i++) randtab[i] = (float)rand()/RAND_MAX;

  // Initialize the simulation state
  reset_simulation();
  microsecond_time(time_program_started);

  glutMainLoop();
}
