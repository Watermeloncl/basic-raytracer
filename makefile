all: ray

ray: ray.c
	gcc -o raytracer ray.c -lm
