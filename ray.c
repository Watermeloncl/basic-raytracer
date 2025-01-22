#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define WINDOW_HEIGHT 2
#define WINDOW_WIDTH 2
#define PIXELS_X 4096
#define PIXELS_Y 4096
#define D 1.2
#define MAX_SPHERES 10
#define MAX_TRIANGLES 12
#define MAX_OBJECT 4194304
#define MAX_VERTICES 10

typedef struct {
	double lightDirX, lightDirY, lightDirZ;
	double red, green, blue;
	double ambRed, ambGreen, ambBlue;
	double backRed, backGreen, backBlue;
} lights_t;

typedef struct {
	double Kd, Ks, Ka;
	double OdRed, OdGreen, OdBlue;
	double OsRed, OsGreen, OsBlue;
	double Kgls;
	double Refl;
	double ambRed, ambGreen, ambBlue;
	double difRed, difGreen, difBlue;
	double specRed, specGreen, specBlue;
} material_t;

typedef struct {
	double center_x, center_y, center_z;
	double radius;
	material_t* mat;
} sphere_t;

typedef struct {
	double vertices[9];
	double nx, ny, nz;
	double rx, ry, rz;
	material_t* mat;
	int inShadow;
} triangle_t;

enum Shape {
	SPHERE,
	TRIANGLE
};

//int times = 0;

void createObjects(unsigned char* fileName, lights_t* lights, sphere_t** spheres, triangle_t** triangles);
void fillObjConsts(lights_t* lights, sphere_t** spheres, triangle_t** triangles);
void fillMaterialConsts(lights_t* lights, material_t* mat);
void fillTriangleConsts(triangle_t* triangle, lights_t* lights);
void assign(unsigned char* pixelArray, int *offset, int fd, unsigned char red, unsigned char green, unsigned char blue);
ssize_t writeToFile(int fd, unsigned char* pixelArray, int numBytes);
void fillPixels(sphere_t** spheres, triangle_t** triangles, lights_t* lights);
void computePixel(double ox, double oy, double oz, double rayDirX, double rayDirY, double rayDirZ, sphere_t** spheres, triangle_t** triangles, enum Shape type, int index, lights_t* lights, unsigned char *red, unsigned char *green, unsigned char *blue);
void calcColor(double px, double py, double pz, double vx, double vy, double vz, double t, lights_t* lights, sphere_t** spheres, triangle_t** triangles, int index, enum Shape currShape, unsigned char *red, unsigned char *green, unsigned char *blue);
int inShadow(double ox, double oy, double oz, sphere_t** spheres, triangle_t** triangles, lights_t* lights, enum Shape type, int index);
void normalizeVector(double *x, double *y, double *z);
double dot(double x1, double y1, double z1, double x2, double y2, double z2);
void cross(double x1, double y1, double z1, double x2, double y2, double z2, double* xf, double* yf, double* zf);
void subVect(double x1, double y1, double z1, double x2, double y2, double z2, double* xf, double* yf, double* zf);
int pointInTriangle(double px, double py, double pz, triangle_t* triangle);
unsigned char convertToPixel(double value);

int main(int argc, char **argv) {

    lights_t* lights = malloc(sizeof(lights_t));
    sphere_t* spheres[MAX_SPHERES];
    triangle_t* triangles[MAX_TRIANGLES];


    for(int i = 0; i < MAX_SPHERES; i++) {
	spheres[i] = NULL;
    }

    for(int i = 0; i < MAX_TRIANGLES; i++) {
	triangles[i] = NULL;
    }

    createObjects("samples/sample1.txt", lights, spheres, triangles);
    fillObjConsts(lights, spheres, triangles);

    fillPixels(spheres, triangles, lights);

    free(lights);
    for(int i = 0; i < MAX_SPHERES; i++) {
	if(spheres[i] != NULL) {
	    free(spheres[i]);
	} else {
	    break;
	}
    }

    for(int i = 0; i < MAX_TRIANGLES; i++) {
	if(triangles[i] != NULL) {
	    free(triangles[i]);
	} else {
	    break;
	}
    }

    exit(0);
}

void createObjects(unsigned char* fileName, lights_t* lights, sphere_t** spheres, triangle_t** triangles) {
    int fd = fileno(fopen(fileName, "r"));

    unsigned char buf[64];
    bzero(buf, 64);

    int nread;
    int offset = 0;
    int values = 0;

    int numSpheres = 0;
    int numTriangles = 0;

    double* valuesToUpdate[32];

    for(;;) {
	nread = read(fd, buf + offset, 1);

	if(buf[offset] == ' ' || buf[offset] == '\n' || nread == 0 || buf[offset] == '\r') {
	    buf[offset] = 0;
	    if(values > 0) {
		if(((buf[0] >= 48) && (buf[0] <= 57)) || (buf[0] == 46) || (buf[0] == 45)) {
		    *(valuesToUpdate[values - 1]) = strtod(buf, NULL);
		    values--;
		}
	    } else if(strcmp(buf, "DirectionToLight") == 0) {
		valuesToUpdate[2] = &(lights->lightDirX);
		valuesToUpdate[1] = &(lights->lightDirY);
		valuesToUpdate[0] = &(lights->lightDirZ);
		values = 3;
	    } else if(strcmp(buf, "LightColor") == 0) {
		valuesToUpdate[2] = &(lights->red);
		valuesToUpdate[1] = &(lights->green);
		valuesToUpdate[0] = &(lights->blue);
		values = 3;
	    } else if(strcmp(buf, "AmbientLight") == 0) {
		valuesToUpdate[2] = &(lights->ambRed);
		valuesToUpdate[1] = &(lights->ambGreen);
		valuesToUpdate[0] = &(lights->ambBlue);
		values = 3;
	    } else if(strcmp(buf, "BackgroundColor") == 0) {
		valuesToUpdate[2] = &(lights->backRed);
		valuesToUpdate[1] = &(lights->backGreen);
		valuesToUpdate[0] = &(lights->backBlue);
		values = 3;
	    } else if((strcmp(buf, "Sphere") == 0) && (numSpheres < MAX_SPHERES - 1)) {
		sphere_t* newSphere = malloc(sizeof(sphere_t));
		spheres[numSpheres] = newSphere;

		material_t* newMaterial = malloc(sizeof(material_t));
		spheres[numSpheres]->mat = newMaterial;

		numSpheres++;

		valuesToUpdate[14] = &(newSphere->center_x);
		valuesToUpdate[13] = &(newSphere->center_y);
		valuesToUpdate[12] = &(newSphere->center_z);
		valuesToUpdate[11] = &(newSphere->radius);
		valuesToUpdate[10] = &(newMaterial->Kd);
		valuesToUpdate[9] = &(newMaterial->Ks);
		valuesToUpdate[8] = &(newMaterial->Ka);
		valuesToUpdate[7] = &(newMaterial->OdRed);
		valuesToUpdate[6] = &(newMaterial->OdGreen);
		valuesToUpdate[5] = &(newMaterial->OdBlue);
		valuesToUpdate[4] = &(newMaterial->OsRed);
		valuesToUpdate[3] = &(newMaterial->OsGreen);
		valuesToUpdate[2] = &(newMaterial->OsBlue);
		valuesToUpdate[1] = &(newMaterial->Kgls);
		valuesToUpdate[0] = &(newMaterial->Refl);

		values = 15;
	    } else if((strcmp(buf, "Triangle") == 0) && (numTriangles < MAX_TRIANGLES - 1)) {
		triangle_t* newTriangle = malloc(sizeof(triangle_t));
		triangles[numTriangles] = newTriangle;

		material_t* newMaterial = malloc(sizeof(material_t));
		triangles[numTriangles]->mat = newMaterial;

		numTriangles++;

		valuesToUpdate[19] = &(newTriangle->vertices[0]);
		valuesToUpdate[18] = &(newTriangle->vertices[1]);
		valuesToUpdate[17] = &(newTriangle->vertices[2]);
		valuesToUpdate[16] = &(newTriangle->vertices[6]);
		valuesToUpdate[15] = &(newTriangle->vertices[7]);
		valuesToUpdate[14] = &(newTriangle->vertices[8]);
		valuesToUpdate[13] = &(newTriangle->vertices[3]);
		valuesToUpdate[12] = &(newTriangle->vertices[4]);
		valuesToUpdate[11] = &(newTriangle->vertices[5]);
		valuesToUpdate[10] = &(newMaterial->Kd);
		valuesToUpdate[9] = &(newMaterial->Ks);
		valuesToUpdate[8] = &(newMaterial->Ka);
		valuesToUpdate[7] = &(newMaterial->OdRed);
		valuesToUpdate[6] = &(newMaterial->OdGreen);
		valuesToUpdate[5] = &(newMaterial->OdBlue);
		valuesToUpdate[4] = &(newMaterial->OsRed);
		valuesToUpdate[3] = &(newMaterial->OsGreen);
		valuesToUpdate[2] = &(newMaterial->OsBlue);
		valuesToUpdate[1] = &(newMaterial->Kgls);
		valuesToUpdate[0] = &(newMaterial->Refl);

		values = 20;
	    }

	    if(nread == 0) {
		break;
	    }

	    bzero(buf, 64);
	    offset = 0;
	} else {
	    offset++;
	}
    }
    close(fd);
}

void fillObjConsts(lights_t* lights, sphere_t** spheres, triangle_t** triangles) {
    normalizeVector(&(lights->lightDirX), &(lights->lightDirY), &(lights->lightDirZ));

    for(int i = 0; i < MAX_SPHERES; i++) {
	if(spheres[i] == NULL) {
	    break;
	}

	fillMaterialConsts(lights, spheres[i]->mat);
    }

    for(int i = 0; i < MAX_TRIANGLES; i++) {
	if(triangles[i] == NULL) {
	    break;
	}

	fillMaterialConsts(lights, triangles[i]->mat);
	fillTriangleConsts(triangles[i], lights);
    }
}

void fillMaterialConsts(lights_t* lights, material_t* mat) {
    mat->ambRed = mat->Ka * lights->ambRed * mat->OdRed;
    mat->ambGreen = mat->Ka * lights->ambGreen * mat->OdGreen;
    mat->ambBlue = mat->Ka * lights->ambBlue * mat->OdBlue;

    mat->difRed = mat->Kd * lights->red * mat->OdRed;
    mat->difGreen = mat->Kd * lights->green * mat->OdGreen;
    mat->difBlue = mat->Kd * lights->blue * mat->OdBlue;

    mat->specRed = mat->Ks * lights->red * mat->OsRed;
    mat->specGreen = mat->Ks * lights->green * mat->OsGreen;
    mat->specBlue = mat->Ks * lights->blue * mat->OsBlue;
}

void fillTriangleConsts(triangle_t* triangle, lights_t* lights) {
    double* vertices = triangle->vertices;

    cross(vertices[0] - vertices[3], vertices[1] - vertices[4],
	  vertices[2] - vertices[5], vertices[6] - vertices[3],
	  vertices[7] - vertices[4], vertices[8] - vertices[5],
	  &(triangle->nx), &(triangle->ny), &(triangle->nz));

    normalizeVector(&(triangle->nx), &(triangle->ny), &(triangle->nz));
    printf("triNormal: (%lf, %lf, %lf)\n", triangle->nx, triangle->ny,
					   triangle->nz);

    double nl = dot(triangle->nx, triangle->ny, triangle->nz, lights->lightDirX, lights->lightDirY, lights->lightDirZ);

    if(nl < 0) {
	triangle->inShadow = 1;
    } else {
	triangle->inShadow = 0;
	triangle->mat->difRed = triangle->mat->difRed * nl;
	triangle->mat->difGreen = triangle->mat->difGreen * nl;
	triangle->mat->difBlue = triangle->mat->difBlue * nl;

	double subR = 2 * dot(triangle->nx, triangle->ny, triangle->nz, lights->lightDirX, lights->lightDirY, lights->lightDirZ);
	triangle->rx = (subR * triangle->nx) - lights->lightDirX;
	triangle->ry = (subR * triangle->ny) - lights->lightDirY;
	triangle->rz = (subR * triangle->nz) - lights->lightDirZ;
    }
}

void assign(unsigned char* pixelArray, int *offset, int fd, unsigned char red, unsigned char green, unsigned char blue) {
    if(*offset >= (MAX_OBJECT - 1)) {
	writeToFile(fd, pixelArray, *offset);
	bzero(pixelArray, MAX_OBJECT);
	*offset = 0;
    }

    *(pixelArray + *offset) = red;
    *(pixelArray + *offset + 1) = green;
    *(pixelArray + *offset + 2) = blue;
    *offset += 3;

}

ssize_t writeToFile(int fd, unsigned char* pixelArray, int numBytes) {

    int nWritten;
    int leftToWrite = numBytes;
    int off = 0;

    for(;;) {
	if(leftToWrite >= MAX_OBJECT) {
	    nWritten = write(fd, pixelArray + off, MAX_OBJECT);
	} else {
	    nWritten = write(fd, pixelArray + off, leftToWrite);
	}

	if(nWritten < 0) {
	    perror("write: ");
	    exit(1);
	} else if((nWritten == 0) && (leftToWrite == 0)) {
	    break;
	} else {
	    leftToWrite -= nWritten;
	    off += nWritten;
	}
    }
}

void fillPixels(sphere_t** spheres, triangle_t** triangles, lights_t* lights) {
    double windowHeight = (double)WINDOW_HEIGHT;
    double windowWidth = (double)WINDOW_WIDTH;
    double pixelsX = (double)PIXELS_X;
    double pixelsY = (double)PIXELS_Y;

    unsigned char color[3];

    double pixelWidth = windowWidth / pixelsX;
    double pixelHeight = windowHeight / pixelsY;

    double startPCX = -(windowWidth / 2) + (0.5 * pixelWidth);
    double startPCY = (windowHeight / 2 ) - (0.5 * pixelHeight);
    double currentPCX = startPCX;
    double currentPCY = startPCY;

    unsigned char pixelArray[MAX_OBJECT];
    bzero(pixelArray, MAX_OBJECT);
    int offset = 0;

    int fd = fileno(fopen("output.ppm", "w"));
    unsigned char buf[32];
    bzero(buf, 32);
    sprintf(buf, "P6\n%d %d\n255\n", PIXELS_X, PIXELS_Y);
    write(fd, buf, strlen(buf));

    for(int y = 0; y < PIXELS_Y; y++) {
	for(int x = 0; x < PIXELS_X; x++) {
	    computePixel(0, 0, (double)D, currentPCX, currentPCY, -(double)D, spheres, triangles, SPHERE, -1, lights, color, color + 1, color + 2);
	    assign(pixelArray, &offset, fd, color[0], color[1], color[2]);

	    currentPCX += pixelWidth;
	}
	currentPCX = startPCX;
	currentPCY -= pixelHeight;
    }

    if(offset > 0) {
	writeToFile(fd, pixelArray, offset);
    }
    close(fd);
}

void computePixel(double ox, double oy, double oz, double rayDirX, double rayDirY, double rayDirZ, sphere_t** spheres, triangle_t** triangles, enum Shape type, int index, lights_t* lights, unsigned char *red, unsigned char *green, unsigned char *blue) {

    normalizeVector(&rayDirX, &rayDirY, &rayDirZ);

    int culled = 0;
    double smallestT = -1;
    double closestIndex = -1;
    enum Shape currShape = SPHERE;
    double px, py, pz;

    for(int i = 0; i < MAX_SPHERES; i++) {
	if(spheres[i] != NULL) {
	    if((type == SPHERE) && (index == i)) {
		continue;
	    }

	    double a = 1;

	    double b = 2 * ((rayDirX * ox) - (rayDirX * spheres[i]->center_x) + (rayDirY * oy) - (rayDirY * spheres[i]->center_y) + (rayDirZ * oz) - (rayDirZ * spheres[i]->center_z));
	    double c = pow(ox, 2) - (2 * ox * spheres[i]->center_x) + pow(spheres[i]->center_x, 2) +
			pow(oy, 2) - (2 * oy * spheres[i]->center_y) + pow(spheres[i]->center_y, 2) +
			pow(oz, 2) - (2 * oz * spheres[i]->center_z) + pow(spheres[i]->center_z, 2) -
			pow(spheres[i]->radius, 2);

	    double disc = pow(b, 2) - (4 * c);
	    double t = -1;
	    if(disc == 0) {
		t = -b/2;
	    } else if(disc > 0) {
		t = (-b - sqrt(pow(b, 2) - (4*c))) / 2;
		if(t < 0) {
		    t = (-b + sqrt(pow(b, 2) - (4*c))) / 2;
		}
	    }

	    if(t > 0) {
		if((smallestT < 0) || (t < smallestT)) {
		    smallestT = t;
		    currShape = SPHERE;
		    closestIndex = i;
		}
	    }

	} else {
	    break;
	}
    }


    for(int i = 0; i < MAX_TRIANGLES; i++) {
	if(triangles[i] != NULL) {
	    if((type == TRIANGLE) && (index == i)) {
		continue;
	    }

	    double t = -1;
	    triangle_t* currTri = triangles[i];

	    double eD = -1 * ((currTri->nx * currTri->vertices[0]) + (currTri->ny * currTri->vertices[1]) + (currTri->nz * currTri->vertices[2]));

	    double nr = dot(currTri->nx, currTri->ny, currTri->nz, rayDirX, rayDirY, rayDirZ);

	    if(nr == 0) {
		continue;
	    }

	    t = (-1 * ((currTri->nx * ox) + (currTri->ny * oy) + (currTri->nz * oz) + eD)) / nr;

	    if(t > 0) {
		if((smallestT < 0) || (t < smallestT)) {
		    double tpx = ox + rayDirX*t;
		    double tpy = oy + rayDirY*t;
		    double tpz = oz + (rayDirZ*t);

		    if(pointInTriangle(tpx, tpy, tpz, currTri)) {
			smallestT = t;
			closestIndex = i;
			currShape = TRIANGLE;
			px = tpx;
			py = tpy;
			pz = tpz;
			if(nr > 0) {
			    culled = 1;
			} else {
			    culled = 0;
			}
		    }
		}
	    }

	} else {
	    break;
	}
    }

    if(smallestT < 0) {
	*red = convertToPixel(lights->backRed);
	*green = convertToPixel(lights->backGreen);
	*blue = convertToPixel(lights->backBlue);
    } else {
	if(currShape == SPHERE) {
	    px = ox + rayDirX*smallestT;
	    py = oy + rayDirY*smallestT;
	    pz = oz + rayDirZ*smallestT;
	} else {
	    if(culled) {
		*red = 0;
		*green = 0;
		*blue = 0;
		return;
	    }
	}

	if(inShadow(px, py, pz, spheres, triangles, lights, currShape, closestIndex)) {
	    *red = 0;
	    *green = 0;
	    *blue = 0;
	} else {
	    calcColor(px, py, pz, rayDirX, rayDirY, rayDirZ, smallestT, lights, spheres, triangles, closestIndex, currShape, red, green, blue);
	}

    }

}

int inShadow(double ox, double oy, double oz, sphere_t** spheres, triangle_t** triangles, lights_t* lights, enum Shape type, int index) {
    double rayDirX = lights->lightDirX;
    double rayDirY = lights->lightDirY;
    double rayDirZ = lights->lightDirZ;

    for(int i = 0; i < MAX_SPHERES; i++) {
	if(spheres[i] != NULL) {
	    if((type == SPHERE) && (index == i)) {
		continue;
	    }

	    double a = 1;

	    double b = 2 * ((rayDirX * ox) - (rayDirX * spheres[i]->center_x) + (rayDirY * oy) - (rayDirY * spheres[i]->center_y) + (rayDirZ * oz) - (rayDirZ * spheres[i]->center_z));
	    double c = pow(ox, 2) - (2 * ox * spheres[i]->center_x) + pow(spheres[i]->center_x, 2) +
			pow(oy, 2) - (2 * oy * spheres[i]->center_y) + pow(spheres[i]->center_y, 2) +
			pow(oz, 2) - (2 * oz * spheres[i]->center_z) + pow(spheres[i]->center_z, 2) -
			pow(spheres[i]->radius, 2);


	    double disc = pow(b, 2) - (4 * c);
	    double t = -1;
	    if(disc == 0) {
		t = -b/2;
	    } else if(disc > 0) {
		t = (-b + sqrt(pow(b, 2) - (4*c))) / 2;
		if(t > 0) {
		    return 1;
		}
	    }

	    if(t > 0) {
		return 1;
	    }

	} else {
	    break;
	}
    }


    for(int i = 0; i < MAX_TRIANGLES; i++) {
	if(triangles[i] != NULL) {
	    if((type == TRIANGLE) && (index == i)) {
		continue;
	    }

	    triangle_t* currTri = triangles[i];

	    double eD = -1 * ((currTri->nx * currTri->vertices[0]) + (currTri->ny * currTri->vertices[1]) + (currTri->nz * currTri->vertices[2]));
	    double nr = dot(currTri->nx, currTri->ny, currTri->nz, rayDirX, rayDirY, rayDirZ);

	    if(nr == 0) {
		continue;
	    }
	    double t = (-1 * ((currTri->nx * ox) + (currTri->ny * oy) + (currTri->nz * oz) + eD)) / nr;

	    if(t > 0) {
		double tpx = ox + rayDirX*t;
		double tpy = oy + rayDirY*t;
		double tpz = oz + (rayDirZ*t);

		if(pointInTriangle(tpx, tpy, tpz, currTri)) {
		    return 1;
		}
	    }

	} else {
	    break;
	}
    }

    return 0;
}



void calcColor(double px, double py, double pz, double vx, double vy, double vz, double t, lights_t* lights, sphere_t** spheres, triangle_t** triangles, int index, enum Shape currShape, unsigned char *red, unsigned char *green, unsigned char *blue) {

    double newRed;
    double newGreen;
    double newBlue;

    material_t* mat;
    double nx;
    double ny;
    double nz;

    if(currShape == TRIANGLE) {
	double vr = dot(vx, vy, vz, triangles[index]->rx, triangles[index]->ry, triangles[index]->rz);
	nx = triangles[index]->nx;
	ny = triangles[index]->ny;
	nz = triangles[index]->nz;
	mat = triangles[index]->mat;

	double vrGls = pow(vr, mat->Kgls);

	if(triangles[index]->inShadow) {
	    newRed = mat->ambRed;
	    newGreen = mat->ambGreen;
	    newBlue = mat->ambBlue;
	} else {
	    newRed = mat->ambRed +
		     mat->difRed +
		     (mat->specRed * vrGls);
	    newGreen = mat->ambGreen +
		       mat->difGreen +
		       (mat->specGreen * vrGls);
	    newBlue = mat->ambBlue +
		      mat->difBlue +
		      (mat->specBlue * vrGls);
	}
    } else {
	mat = spheres[index]->mat;

	nx = (px - spheres[index]->center_x) / spheres[index]->radius;
	ny = (py - spheres[index]->center_y) / spheres[index]->radius;
	nz = (pz - spheres[index]->center_z) / spheres[index]->radius;

	double nl = dot(nx, ny, nz, lights->lightDirX, lights->lightDirY, lights->lightDirZ);

	if(nl <= 0) {
	    newRed = mat->ambRed;
	    newGreen = mat->ambGreen;
	    newBlue = mat->ambBlue;
	} else {
	    double rx = (2 * nl * nx) - lights->lightDirX;
	    double ry = (2 * nl * ny) - lights->lightDirY;
	    double rz = (2 * nl * nz) - lights->lightDirZ;
	    double vr = dot(vx, vy, vz, rx, ry, rz);

	    newRed = mat->ambRed +
		 (mat->difRed * nl) +
		 (mat->specRed * pow(vr, mat->Kgls));
	    newGreen = mat->ambGreen +
		   (mat->difGreen * nl) +
		   (mat->specGreen * pow(vr, mat->Kgls));

	    newBlue = mat->ambBlue +
		   (mat->difBlue * nl) +
		   (mat->specBlue * pow(vr, mat->Kgls));
	}
    }

    if(mat->Refl > 0) {
	unsigned char tempRed;
	unsigned char tempGreen;
	unsigned char tempBlue;

	double newDirX;
	double newDirY;
	double newDirZ;

	double subNewDir = 2 * dot(vx, vy, vz, nx, ny, nz);

	newDirX = vx - (nx * subNewDir);
	newDirY = vy - (ny * subNewDir);
	newDirZ = vz - (nz * subNewDir);
	normalizeVector(&newDirX, &newDirY, &newDirZ);

	computePixel(px, py, pz, newDirX, newDirY, newDirZ, spheres, triangles, currShape, index, lights, &tempRed, &tempGreen, &tempBlue);

/*	//this is what it should be!
	*red = ((1 - mat->Refl) * convertToPixel(newRed)) + (mat->Refl * tempRed);
	*green = ((1 - mat->Refl) * convertToPixel(newGreen)) + (mat->Refl * tempGreen);
	*blue = ((1 - mat->Refl) * convertToPixel(newBlue)) + (mat->Refl * tempBlue);
*/
	//this is what's in the slides, and also in the example images, so my pics 4/5 follow suit
	*red = convertToPixel(newRed) + (mat->Refl * tempRed);
	*green = convertToPixel(newGreen) + (mat->Refl * tempGreen);
	*blue = convertToPixel(newBlue) + (mat->Refl * tempBlue);

    } else {
	*red = convertToPixel(newRed);
	*green = convertToPixel(newGreen);
	*blue = convertToPixel(newBlue);
    }
}

void normalizeVector(double *x, double *y, double *z) {
    double mag = sqrt((pow(*x, 2.0) + pow(*y, 2.0) + pow(*z, 2.0)));
    *x = *x / mag;
    *y = *y / mag;
    *z = *z / mag;
}

double dot(double x1, double y1, double z1, double x2, double y2, double z2) {
    return (x1 * x2) + (y1 * y2) + (z1 * z2);
}

void cross(double x1, double y1, double z1, double x2, double y2, double z2, double* xf, double* yf, double* zf) {
    *xf = (y1 * z2) - (z1 * y2);
    *yf = (z1 * x2) - (x1 * z2);
    *zf = (x1 * y2) - (y1 * x2);
}

void subVect(double x1, double y1, double z1, double x2, double y2, double z2, double* xf, double* yf, double* zf) {
    *xf = x1 - x2;
    *yf = y1 - y2;
    *zf = z1 - z2;
}

int pointInTriangle(double px, double py, double pz, triangle_t* triangle) {
    double* vertices = triangle->vertices;
    double result;
    double xf, yf, zf;


    cross(px - vertices[3], py - vertices[4], pz - vertices[5],
	  vertices[6] - vertices[3], vertices[7] - vertices[4],
	  vertices[8] - vertices[5], &xf, &yf, &zf);

    result = dot(triangle->nx, triangle->ny, triangle->nz, xf, yf, zf);

    if(result < 0) {
	return 0;
    }

    cross(px - vertices[6], py - vertices[7], pz - vertices[8],
	  vertices[0] - vertices[6], vertices[1] - vertices[7],
	  vertices[2] - vertices[8], &xf, &yf, &zf);

    result = dot(triangle->nx, triangle->ny, triangle->nz, xf, yf, zf);

    if(result < 0) {
	return 0;
    }

    cross(px - vertices[0], py - vertices[1], pz - vertices[2],
	  vertices[3] - vertices[0], vertices[4] - vertices[1],
	  vertices[5] - vertices[2], &xf, &yf, &zf);

    result = dot(triangle->nx, triangle->ny, triangle->nz, xf, yf, zf);

    if(result < 0) {
	return 0;
    }


    return 1;
}

unsigned char convertToPixel(double value) {
    if(value > 1) {
	value = 1;
    } else if(value < 0) {
	value = 0;
    }

    return (unsigned char)round(255*value);
}
