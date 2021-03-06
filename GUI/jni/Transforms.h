#ifndef TRANSFORMS_H
#define TRANSFORMS_H
#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <math.h>

typedef struct Transform {
		int points;
		float* sine;
		float* cosine;
		float* real;
		float* imaginary;
		void (*doTransform)(struct Transform* transform, float* input);
} Transform;

Transform* newTransform(int type, int points);
void transformMagnitude(Transform* transform, float* output);
void destroyTransform(Transform** transform);

#endif
