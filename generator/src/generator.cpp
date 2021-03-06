// generator.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <cstring>
#include <sstream>

using namespace std;

FILE *ficheiro;

void printErrorArgs(char* func) {

	printf("[Shape **%s**]: Invalid Arguments! \n",func);
}


void plane(float x) {
		
		fprintf(ficheiro, "%d\n", 6);

		if ((x = fabs(x)) == 0) x = 1;

		fprintf(ficheiro, "%f %f %f\n", -x / 2, 0.0, x / 2);
		fprintf(ficheiro, "%f %f %f\n",  x / 2, 0.0, x / 2);
		fprintf(ficheiro, "%f %f %f\n",  x / 2, 0.0, -x / 2);

		fprintf(ficheiro, "%f %f %f\n", x / 2, 0.0, -x / 2);
		fprintf(ficheiro, "%f %f %f\n", -x / 2, 0.0, -x / 2);
		fprintf(ficheiro, "%f %f %f\n", -x / 2, 0.0, x / 2);


}


void box(float x, float z, float y, int divisions = 0) {


	if ((divisions = abs(divisions)) == 0) divisions = 1;

	fprintf(ficheiro, "%ld\n", divisions*divisions*6*6);

	float xaux = -x / 2;
	float zaux = z / 2;
	float yaux = y / 2;

	//frente


	for (int j = 0; j < divisions; j++) {

		for (int i = 0; i < divisions; i++) {

			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, z / 2);
			fprintf(ficheiro, "%f %f %f\n", xaux, yaux - (y / divisions), z / 2);
			fprintf(ficheiro, "%f %f %f\n", xaux + (x / divisions), yaux - (y / divisions), z / 2);

			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, z / 2);
			fprintf(ficheiro, "%f %f %f\n", xaux + (x / divisions), yaux - (y / divisions), z / 2);
			fprintf(ficheiro, "%f %f %f\n", xaux + (x / divisions), yaux, z / 2);

			xaux += (x / divisions);
		}

		xaux = -x / 2;
		yaux = yaux - (y / divisions);
	}



	yaux = y / 2;
	xaux = -x / 2;
	zaux = z / 2;

	//atras

	for (int j = 0; j < divisions; j++) {

		for (int i = 0; i < divisions; i++) {

			fprintf(ficheiro, "%f %f %f\n", xaux + (x / divisions), yaux, -z / 2);
			fprintf(ficheiro, "%f %f %f\n", xaux + (x / divisions), yaux - (y / divisions), -z / 2);
			fprintf(ficheiro, "%f %f %f\n", xaux, yaux - (y / divisions), -z / 2);

			fprintf(ficheiro, "%f %f %f\n", xaux, yaux - (y / divisions), -z / 2);
			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, -z / 2);
			fprintf(ficheiro, "%f %f %f\n", xaux + (x / divisions), yaux, -z / 2);

			xaux += (x / divisions);
		}

		xaux = -x / 2;
		yaux = yaux - (y / divisions);
	}


	xaux = x / 2;
	zaux = z / 2;
	yaux = y / 2;

	//dir

	for (int j = 0; j < divisions; j++) {

		for (int i = 0; i < divisions; i++) {

			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, zaux);
			fprintf(ficheiro, "%f %f %f\n", xaux, yaux - (y / divisions), zaux);
			fprintf(ficheiro, "%f %f %f\n", xaux, yaux - (y / divisions), zaux - (z / divisions));

			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, zaux);
			fprintf(ficheiro, "%f %f %f\n", xaux, yaux - (y / divisions), zaux - (z / divisions));
			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, zaux - (z / divisions));

			zaux -= (z / divisions);
		}

		zaux = z / 2;
		yaux = yaux - (y / divisions);
	}


	xaux = -x / 2;
	zaux = z / 2;
	yaux = y / 2;

	//esq

	for (int j = 0; j < divisions; j++) {

		for (int i = 0; i < divisions; i++) {


			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, zaux - (z / divisions));
			fprintf(ficheiro, "%f %f %f\n", xaux, yaux - (y / divisions), zaux);
			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, zaux);

			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, zaux - (z / divisions));
			fprintf(ficheiro, "%f %f %f\n", xaux, yaux - (y / divisions), zaux - (z / divisions));
			fprintf(ficheiro, "%f %f %f\n", xaux, yaux - (y / divisions), zaux);

			zaux -= (z / divisions);
		}

		zaux = z / 2;
		yaux = yaux - (y / divisions);
	}



	xaux = -x / 2;
	zaux = -z / 2;
	yaux = y / 2;

	//topo

	for (int j = 0; j < divisions; j++) {

		for (int i = 0; i < divisions; i++) {

			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, zaux);
			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, zaux + (z / divisions));
			fprintf(ficheiro, "%f %f %f\n", xaux + (x / divisions), yaux, zaux + (z / divisions));

			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, zaux);
			fprintf(ficheiro, "%f %f %f\n", xaux + (x / divisions), yaux, zaux + (z / divisions));
			fprintf(ficheiro, "%f %f %f\n", xaux + (x / divisions), yaux, zaux);

			xaux += (x / divisions);
		}

		xaux = -x / 2;
		zaux = zaux + (y / divisions);
	}


	xaux = -x / 2;
	zaux = -z / 2;
	yaux = -y / 2;

	//base

	for (int j = 0; j < divisions; j++) {

		for (int i = 0; i < divisions; i++) {

			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, zaux);
			fprintf(ficheiro, "%f %f %f\n", xaux + (x / divisions), yaux, zaux + (z / divisions));
			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, zaux + (z / divisions));

			fprintf(ficheiro, "%f %f %f\n", xaux, yaux, zaux);
			fprintf(ficheiro, "%f %f %f\n", xaux + (x / divisions), yaux, zaux);
			fprintf(ficheiro, "%f %f %f\n", xaux + (x / divisions), yaux, zaux + (z / divisions));

			xaux += (x / divisions);
		}

		xaux = -x / 2;
		zaux = zaux + (y / divisions);
	}
}


void cone(float radius, float height, int slices, int stacks) {

	//glPushMatrix();
	//glTranslatef(0, height / 2 - radius, radius);
	//glRotatef(90, 1, 0, 0);


	fprintf(ficheiro, "%ld\n", (6*slices)+(6*slices*(stacks-1)));


	float angulo = (2 * M_PI) / slices;
	float baseaux = -height / 2;
	float raioaux = radius;
	float raioaux2 = radius - (radius / stacks);

	//DESENHA BASE


	for (int i = 0; i<slices; i++) {
		fprintf(ficheiro, "%f %f %f\n", radius*sin(angulo*i), baseaux, radius*cos(angulo*i));
		fprintf(ficheiro, "%f %f %f\n", 0.0f, baseaux, 0.0f);
		fprintf(ficheiro, "%f %f %f\n", radius*sin(angulo*(i + 1)), baseaux, radius*cos(angulo*(i + 1)));
	}



	//DESENHA STACKS
	for (int j = 0; j < stacks-1; j++) {

		for (int i = 0; i < slices; i++) {

			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*i)), baseaux, raioaux*(cos(angulo*i)));
			fprintf(ficheiro, "%f %f %f\n", raioaux2*(sin(angulo*(i + 1))), baseaux + (height / stacks), raioaux2*(cos(angulo*(i + 1))));
			fprintf(ficheiro, "%f %f %f\n", raioaux2*(sin(angulo*i)), baseaux + (height / stacks), raioaux2*(cos(angulo*i)));

			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*i)), baseaux, raioaux*(cos(angulo*i)));
			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*(i + 1))), baseaux, raioaux*(cos(angulo*(i + 1))));
			fprintf(ficheiro, "%f %f %f\n", raioaux2*(sin(angulo*(i + 1))), baseaux + (height / stacks), raioaux2*(cos(angulo*(i + 1))));


		}

		baseaux += height / stacks;
		raioaux = raioaux2;
		raioaux2 -= (radius / stacks);
	}
	
	for (int i = 0; i < slices; i++) {
		fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*i)), baseaux, raioaux*(cos(angulo*i)));
		fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*(i + 1))), baseaux, raioaux*(cos(angulo*(i + 1))));
		fprintf(ficheiro, "%f %f %f\n", raioaux2*(sin(angulo*(i + 1))), baseaux + (height / stacks), raioaux2*(cos(angulo*(i + 1))));
	}


	//glPopMatrix();
}


void sphere(float radius, int slices, int stacks) {


	float alfa = M_PI / stacks;
	float alfaux = M_PI / stacks;
	float altura = 0;
	float altura2 = sin(alfa) * radius;
	float angulo = (2 * M_PI) / slices;
	float raioaux = radius;
	float raioaux2 = sqrtf((radius*radius) - (altura2 * altura2));

	long int num_vertices = ((stacks - 2)*(6 * slices)) + (6 * slices);


	int vertex = 0;

	float *n, *t;
	n = (float *)malloc(sizeof(float) * num_vertices * 3);
	t = (float *)malloc(sizeof(float) * num_vertices * 2);

	fprintf(ficheiro, "%ld\n", num_vertices);

	if (stacks % 2 == 1) {


		alfa = alfa / 2;
		altura = sin(alfa) * radius;
		altura2 = sin(alfa + alfaux) * radius;
		raioaux = sqrtf((radius*radius) - (altura * altura));;
		raioaux2 = sqrtf((radius*radius) - (altura2 * altura2));


		for (int i = 0; i < slices; i++) {


			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*i)), -altura, raioaux*((cos(angulo*i))));
			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*(i + 1))), -altura, raioaux*((cos(angulo*(i + 1)))));
			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*(i + 1))), altura, raioaux*((cos(angulo*(i + 1)))));

			n[vertex * 3 + 0] = raioaux * sin(angulo*i) / radius;
			n[vertex * 3 + 1] = -altura / radius;
			n[vertex * 3 + 2] = raioaux * cos(angulo*i) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			n[vertex * 3 + 0] = raioaux * sin(angulo*(i + 1)) / radius;
			n[vertex * 3 + 1] = -altura / radius;
			n[vertex * 3 + 2] = raioaux * cos(angulo*(i + 1)) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			n[vertex * 3 + 0] = raioaux * sin(angulo*(i + 1)) / radius;
			n[vertex * 3 + 1] = altura / radius;
			n[vertex * 3 + 2] = raioaux * cos(angulo*(i + 1)) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*i)), -altura, raioaux*((cos(angulo*i))));
			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*(i + 1))), altura, raioaux*((cos(angulo*(i + 1)))));
			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*(i))), altura, raioaux*((cos(angulo*(i)))));

			n[vertex * 3 + 0] = raioaux * sin(angulo*i) / radius;
			n[vertex * 3 + 1] = -altura / radius;
			n[vertex * 3 + 2] = raioaux * cos(angulo*i) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			n[vertex * 3 + 0] = raioaux * sin(angulo*(i + 1)) / radius;
			n[vertex * 3 + 1] = altura / radius;
			n[vertex * 3 + 2] = raioaux * cos(angulo*(i + 1)) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			n[vertex * 3 + 0] = raioaux * sin(angulo*(i)) / radius;
			n[vertex * 3 + 1] = altura / radius;
			n[vertex * 3 + 2] = raioaux * cos(angulo*(i)) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;


		}


		alfa += alfaux;
		stacks--;
		
	}


	//DESENHA STACKS
	for (int j = 0; j < (stacks / 2)-1; j++) {

		for (int i = 0; i < slices; i++) {

			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*i)), altura, raioaux*((cos(angulo*i))));
			fprintf(ficheiro, "%f %f %f\n", raioaux2*(sin(angulo*(i + 1))), altura2, raioaux2*((cos(angulo*(i + 1)))));
			fprintf(ficheiro, "%f %f %f\n", raioaux2*(sin(angulo*(i))), altura2, raioaux2*((cos(angulo*(i)))));
			
			n[vertex * 3 + 0] = raioaux * sin(angulo*i) / radius;
			n[vertex * 3 + 1] = altura / radius;
			n[vertex * 3 + 2] = raioaux * cos(angulo*i) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex * 3 + 1]) / M_PI);
			vertex++;

			n[vertex * 3 + 0] = raioaux2 * sin(angulo*(i + 1)) / radius;
			n[vertex * 3 + 1] = altura2 / radius;
			n[vertex * 3 + 2] = raioaux2 * cos(angulo*(i + 1)) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			n[vertex * 3 + 0] = raioaux2 * sin(angulo*(i)) / radius;
			n[vertex * 3 + 1] = altura2 / radius;
			n[vertex * 3 + 2] = raioaux2 * cos(angulo*(i)) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*i)), altura, raioaux*((cos(angulo*i))));
			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*(i + 1))), altura, raioaux*((cos(angulo*(i + 1)))));
			fprintf(ficheiro, "%f %f %f\n", raioaux2*(sin(angulo*(i + 1))), altura2, raioaux2*((cos(angulo*(i + 1)))));
			
			n[vertex * 3 + 0] = raioaux * sin(angulo*i) / radius;
			n[vertex * 3 + 1] = altura / radius;
			n[vertex * 3 + 2] = raioaux * cos(angulo*i) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			n[vertex * 3 + 0] = raioaux * sin(angulo*(i + 1)) / radius;
			n[vertex * 3 + 1] = altura / radius;
			n[vertex * 3 + 2] = raioaux * cos(angulo*(i + 1)) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			n[vertex * 3 + 0] = raioaux2 * sin(angulo*(i+1)) / radius;
			n[vertex * 3 + 1] = altura2 / radius;
			n[vertex * 3 + 2] = raioaux2 * cos(angulo*(i+1)) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;
		}

		for (int i = 0; i < slices; i++) {

			fprintf(ficheiro, "%f %f %f\n", raioaux2*(sin(angulo*i)), -altura2, raioaux2*((cos(angulo*i))));
			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*(i + 1))), -altura, raioaux*((cos(angulo*(i + 1)))));
			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*(i))), -altura, raioaux*((cos(angulo*(i)))));

			n[vertex * 3 + 0] = raioaux2 * sin(angulo*i) / radius;
			n[vertex * 3 + 1] = -altura2 / radius;
			n[vertex * 3 + 2] = raioaux2 * cos(angulo*i) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			n[vertex * 3 + 0] = raioaux * sin(angulo*(i + 1)) / radius;
			n[vertex * 3 + 1] = -altura / radius;
			n[vertex * 3 + 2] = raioaux * cos(angulo*(i + 1)) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			n[vertex * 3 + 0] = raioaux * sin(angulo*(i)) / radius;
			n[vertex * 3 + 1] = -altura / radius;
			n[vertex * 3 + 2] = raioaux * cos(angulo*(i)) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			fprintf(ficheiro, "%f %f %f\n", raioaux2*(sin(angulo*i)), -altura2, raioaux2*((cos(angulo*i))));
			fprintf(ficheiro, "%f %f %f\n", raioaux2*(sin(angulo*(i + 1))), -altura2, raioaux2*((cos(angulo*(i + 1)))));
			fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*(i + 1))), -altura, raioaux*((cos(angulo*(i + 1)))));

			n[vertex * 3 + 0] = raioaux2 * sin(angulo*i) / radius;
			n[vertex * 3 + 1] = -altura2 / radius;
			n[vertex * 3 + 2] = raioaux2 * cos(angulo*i) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			n[vertex * 3 + 0] = raioaux2 * sin(angulo*(i + 1)) / radius;
			n[vertex * 3 + 1] = -altura2 / radius;
			n[vertex * 3 + 2] = raioaux2 * cos(angulo*(i + 1)) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

			n[vertex * 3 + 0] = raioaux * sin(angulo*(i+1)) / radius;
			n[vertex * 3 + 1] = -altura / radius;
			n[vertex * 3 + 2] = raioaux * cos(angulo*(i+1)) / radius;
			t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
			t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
			vertex++;

		}

		alfa += alfaux;
		altura = altura2;
		altura2 = sin(alfa) * radius;
		raioaux = raioaux2;
		raioaux2 = sqrtf((radius*radius) - (altura2 * altura2));
	}

	for (int i = 0; i < slices; i++) {

		fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*i)), altura, raioaux*((cos(angulo*i))));
		fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*(i + 1))), altura, raioaux*((cos(angulo*(i + 1)))));
		fprintf(ficheiro, "%f %f %f\n", raioaux2*(sin(angulo*(i + 1))), altura2, raioaux2*((cos(angulo*(i + 1)))));

		n[vertex * 3 + 0] = raioaux * sin(angulo*i) / radius;
		n[vertex * 3 + 1] = altura / radius;
		n[vertex * 3 + 2] = raioaux * cos(angulo*i) / radius;
		t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
		t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
		vertex++;

		n[vertex * 3 + 0] = raioaux * sin(angulo*(i + 1)) / radius;
		n[vertex * 3 + 1] = altura / radius;
		n[vertex * 3 + 2] = raioaux * cos(angulo*(i + 1)) / radius;
		t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
		t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
		vertex++;

		n[vertex * 3 + 0] = raioaux2 * sin(angulo*(i+1)) / radius;
		n[vertex * 3 + 1] = altura2 / radius;
		n[vertex * 3 + 2] = raioaux2 * cos(angulo*(i+1)) / radius;
		t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
		t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
		vertex++;

		fprintf(ficheiro, "%f %f %f\n", raioaux2*(sin(angulo*i)), -altura2, raioaux2*((cos(angulo*i))));
		fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*(i + 1))), -altura, raioaux*((cos(angulo*(i + 1)))));
		fprintf(ficheiro, "%f %f %f\n", raioaux*(sin(angulo*(i))), -altura, raioaux*((cos(angulo*(i)))));

		n[vertex * 3 + 0] = raioaux2 * sin(angulo*i) / radius;
		n[vertex * 3 + 1] = -altura2 / radius;
		n[vertex * 3 + 2] = raioaux2 * cos(angulo*i) / radius;
		t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
		t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
		vertex++;

		n[vertex * 3 + 0] = raioaux * sin(angulo*(i + 1)) / radius;
		n[vertex * 3 + 1] = -altura / radius;
		n[vertex * 3 + 2] = raioaux * cos(angulo*(i + 1)) / radius;
		t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
		t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
		vertex++;

		n[vertex * 3 + 0] = raioaux * sin(angulo*(i))/radius;
		n[vertex * 3 + 1] = -altura / radius;
		n[vertex * 3 + 2] = raioaux * cos(angulo*(i)) / radius;
		t[vertex * 2 + 0] = atan2(n[vertex * 3 + 2]/radius, n[vertex * 3 + 0]/radius) / (2 * M_PI) + 0.5;
		t[vertex * 2 + 1] = 0.5-(asin(n[vertex*3+1])/M_PI);
		vertex++;
	}
	

	for (int i = 0; i < num_vertices; i++) {
		fprintf(ficheiro, "%f %f %f\n", n[3*i], n[i*3+1], n[i*3+2]);
	}

	for (int i = 0; i < num_vertices; i++) {
		fprintf(ficheiro, "%f %f\n", t[2*i], t[i*2+1]);
	}

	free(t);
	free(n);

}









#define POINT_COUNT 16
// Points that make up the loop for catmull-rom interpolation
/*
void buildRotMatrix(float *x, float *y, float *z, float *m) {

	m[0] = x[0]; m[1] = x[1]; m[2] = x[2]; m[3] = 0;
	m[4] = y[0]; m[5] = y[1]; m[6] = y[2]; m[7] = 0;
	m[8] = z[0]; m[9] = z[1]; m[10] = z[2]; m[11] = 0;
	m[12] = 0; m[13] = 0; m[14] = 0; m[15] = 1;
}*/


void cross(float *a, float *b, float *res) {

	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
}


void normalize(float *a) {

	float l = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	a[0] = a[0] / l;
	a[1] = a[1] / l;
	a[2] = a[2] / l;
}


float length(float *v) {

	float res = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	return res;

}

void multMatrixVector(float *m, float *v, float *res) {

	for (int j = 0; j < 4; ++j) {
		res[j] = 0;
		for (int k = 0; k < 4; ++k) {
			res[j] += v[k] * m[j * 4 + k];
		}
	}

}


float p[4][3] = {};



void getCatmullRomPoint(float t, float *p0, float *p1, float *p2, float *p3, float *pos, float *deriv) {

	// catmull-rom matrix
	float m[4][4] = { { -0.5f,  1.5f, -1.5f,  0.5f },
	{ 1.0f, -2.5f,  2.0f, -0.5f },
	{ -0.5f,  0.0f,  0.5f,  0.0f },
	{ 0.0f,  1.0f,  0.0f,  0.0f } };

	// Compute A = M * P
	float a[4][4];

	float px[4] = { p0[0], p1[0], p2[0], p3[0] };
	float py[4] = { p0[1], p1[1], p2[1], p3[1] };
	float pz[4] = { p0[2], p1[2], p2[2], p3[2] };
	float pw[4] = { p0[3], p1[3], p2[3], p3[3] };

	multMatrixVector((float *)m, px, a[0]);
	multMatrixVector((float *)m, py, a[1]);
	multMatrixVector((float *)m, pz, a[2]);
	multMatrixVector((float *)m, pw, a[3]);

	float tv[4] = { powf(t, 3.0), powf(t, 2.0), t, 1.0 };
	float tdv[4] = { 3.0 * powf(t, 2.0), 2.0f * t, 1.0, 0.0 };

	// Compute pos = T * A
	multMatrixVector((float *)a, tv, pos);
	// Compute deriv = T' * A
	multMatrixVector((float *)a, tdv, deriv);
}


// given  global t, returns the point in the curve
void getGlobalCatmullRomPoint(float gt, float *pos, float *deriv) {



	float t = gt * POINT_COUNT; // this is the real global t
	int index = floor(t);  // which segment
	t = t - index; // where within  the segment

				   // indices store the points
	int indices[4];
	indices[0] = (index + POINT_COUNT - 1) % POINT_COUNT;
	indices[1] = (indices[0] + 1) % POINT_COUNT;
	indices[2] = (indices[1] + 1) % POINT_COUNT;
	indices[3] = (indices[2] + 1) % POINT_COUNT;

	getCatmullRomPoint(t, p[indices[0]], p[indices[1]], p[indices[2]], p[indices[3]], pos, deriv);
}






void read_patch(FILE *f_patch, FILE *f, int tess) {


	int n_patches, n_vertices, k;


	fscanf_s(f_patch, "%d\n", &n_patches);
	
	int** patches = new int*[n_patches];
	for (int i = 0; i < n_patches; i++)
		patches[i] = new int[16];



	for (int i = 0; i < n_patches; i++) {
		int j;
		for (j = 0; j < 15; j++) {

			fscanf_s(f_patch, "%d, ", &k);
			patches[i][j] = k;
		}

		fscanf_s(f_patch, "%d\n", &k);
		patches[i][j] = k;
	}



	fscanf_s(f_patch, "%d\n", &n_vertices);

	float** vertices = new float*[n_vertices];
	for (int i = 0; i < n_vertices; i++)
		vertices[i] = new float[3];




	float xx, yy, zz;
	int v = 0;

	while (fscanf_s(f_patch, " %f, %f, %f\n", &xx, &yy, &zz) != EOF) {
		vertices[v][0] = xx ;
		vertices[v][1] = yy;
		vertices[v++][2] = zz;
	}



	int count = 1;

	float t = 0;

	for (int a = 0; a < n_patches; a++) {
		
		float pg[16][3];


		for (int i = 0; i < 16; i++) {
			pg[i][0] = vertices[patches[a][i]][0];
			pg[i][1] = vertices[patches[a][i]][1];
			pg[i][2] = vertices[patches[a][i]][2];

		}


		for (int q = 0; q < 4; q++) {


			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 3; j++) {
					p[i][j] = pg[4 * q + i][j];	
				}
			}

			while (t <= 1) {

				static float up[4] = { 0, 1, 0, 0 };

				float pos[4], deriv[4];
				float m[4 * 4];
				float y[4], z[4];

				getGlobalCatmullRomPoint(t, pos, deriv);

				cross(deriv, up, z);
				cross(z, deriv, y);

				normalize(deriv);
				normalize(y);
				normalize(z);

				fprintf(f, "%f %f %f\n", pos[0], pos[1], pos[2]);

				up[0] = y[0]; up[1] = y[1]; up[2] = y[2]; up[3] = y[3];

				t += (1.0 / tess);
			}


			t = 0;

		}

	}

}



void generate(char* function, char** args, int argc) {


	if (strcmp(function, "bezier") == 0) {

		string s = args[0];

		string function3d = (string)s.substr(0, s.find(".")) + ".3d";

		fopen_s(&ficheiro, function3d.c_str(), "w");

		FILE *patch;
		if (fopen_s(&patch, args[0], "r") == 0) {

			read_patch(patch, ficheiro, atoi(args[1]));
		}
		else
			printf("Erro ao carregar Bezier control points.");

		fclose(ficheiro);
		fclose(patch);
		return;
	}


	string function3d = (string)function + ".3d";

	fopen_s(&ficheiro, function3d.c_str(), "w");



	if (strcmp(function, "plane") == 0) {

		float x;
		if (argc == 0) x = 0.0;
		else x = atof(args[0]);
		plane(x);
	}

	if (strcmp(function, "box") == 0) {

		int divs;
		if (argc == 0) divs = 1;
		else divs = atoi(args[3]);

		box(atof(args[0]), atof(args[1]), atof(args[2]), divs);
	}

	if (strcmp(function, "cone") == 0) {

		cone(atof(args[0]), atof(args[1]), atoi(args[2]), atoi(args[3]));
	}

	if (strcmp(function, "sphere") == 0) {

		sphere(atof(args[0]), atoi(args[1]), atoi(args[2]));
	}

	fclose(ficheiro);
}






int main(int argc, char** argv)  {

	if (argc < 2) {
		printf("Main: Lack of arguments\n");
		_exit(0);
	}


	if (strcmp(argv[1], "bezier") == 0) {
		if (argc != 4) {
			printErrorArgs(argv[1]);
			_exit(0);
		}

		printf("bezier patch\n");

		generate(argv[1], argv + 2, argc - 2);

		return 0;
	}



	if (strcmp(argv[1], "plane") == 0) {
		if (argc != 2 && argc !=3) {
			printErrorArgs(argv[1]);
			_exit(0);
		}

		printf("plane\n");

		generate(argv[1], argv + 2, argc - 2);

		return 0;
	}
	
	if (strcmp(argv[1], "box") == 0) {

		if (argc != 5 && argc != 6) {
			printErrorArgs(argv[1]);
			_exit(0);
		}
		printf("box\n");

		generate(argv[1], argv + 2, argc - 2);

		return 0;
		
	}

	if (strcmp(argv[1], "cone") == 0) {

		if (argc != 6) {
			printErrorArgs(argv[1]);
			_exit(0);
		}
		printf("cone\n");

		generate(argv[1], argv + 2, argc - 2);

		return 0;
	}

	if (strcmp(argv[1], "sphere") == 0) {

		if (argc != 5) {
			printErrorArgs(argv[1]);
			_exit(0);
		}
		printf("sphere\n");

		generate(argv[1], argv + 2, argc - 2);

		return 0;
	}


	printf("Main: Shape not found!\n");

	return -1;

}

