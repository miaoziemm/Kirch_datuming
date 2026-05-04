#include "../include/se_alloc.h"

/* allocate a 1-d array */
void *alloc1 (size_t n1, size_t size)
{
    void *p;

    if ((p=malloc(n1*size))==NULL)
	return NULL;
    return p;
}

/* re-allocate a 1-d array */
void *realloc1(void *v, size_t n1, size_t size)
{
    void *p;

    if ((p=realloc(v,n1*size))==NULL)
	return NULL;
    return p;
}

/* free a 1-d array */
void free1 (void *p)
{
    free(p);
}

/* allocate a 2-d array */
void **alloc2 (size_t n1, size_t n2, size_t size)
{
    size_t i2;
    void **p;

    if ((p=(void**)malloc(n2*sizeof(void*)))==NULL)
	return NULL;
    if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
        free(p);
        return NULL;
    }
    for (i2=0; i2<n2; i2++)
	p[i2] = (char*)p[0]+size*n1*i2;
    return p;
}

/* free a 2-d array */
void free2 (void **p)
{
    free(p[0]);
    free(p);
}

/* allocate a 3-d array */
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size)
{
    size_t i3,i2;
    void ***p;

    if ((p=(void***)malloc(n3*sizeof(void**)))==NULL)
	return NULL;
    if ((p[0]=(void**)malloc(n3*n2*sizeof(void*)))==NULL) {
        free(p);
        return NULL;
    }
    if ((p[0][0]=(void*)malloc(n3*n2*n1*size))==NULL) {
        free(p[0]);
        free(p);
        return NULL;
    }

    for (i3=0; i3<n3; i3++) {
        p[i3] = p[0]+n2*i3;
        for (i2=0; i2<n2; i2++)
	    p[i3][i2] = (char*)p[0][0]+size*n1*(i2+n2*i3);
    }
    return p;
}

/* free a 3-d array */
void free3 (void ***p)
{
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/* allocate a 4-d array */
void ****alloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size)
{
    size_t i4,i3,i2;
    void ****p;

    if ((p=(void****)malloc(n4*sizeof(void***)))==NULL)
	return NULL;
    if ((p[0]=(void***)malloc(n4*n3*sizeof(void**)))==NULL) {
        free(p);
        return NULL;
    }
    if ((p[0][0]=(void**)malloc(n4*n3*n2*sizeof(void*)))==NULL) {
        free(p[0]);
        free(p);
        return NULL;
    }
    if ((p[0][0][0]=(void*)malloc(n4*n3*n2*n1*size))==NULL) {
        free(p[0][0]);
        free(p[0]);
        free(p);
        return NULL;
    }
    for (i4=0; i4<n4; i4++) {
        p[i4] = p[0]+i4*n3;
        for (i3=0; i3<n3; i3++) {
            p[i4][i3] = p[0][0]+n2*(i3+n3*i4);
            for (i2=0; i2<n2; i2++)
		p[i4][i3][i2] = (char*)p[0][0][0]+
		    size*n1*(i2+n2*(i3+n3*i4));
        }
    }
    return p;
}

/* free a 4-d array */
void free4 (void ****p)
{
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/* allocate a 5-d array */
void *****alloc5 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size)
{
    size_t i5,i4,i3,i2;
    void *****p;

    if ((p=(void*****)malloc(n5*sizeof(void****)))==NULL)
	return NULL;
    if ((p[0]=(void****)malloc(n5*n4*sizeof(void***)))==NULL) {
        free(p);
        return NULL;
    }
    if ((p[0][0]=(void***)malloc(n5*n4*n3*sizeof(void**)))==NULL) {
        free(p[0]);
        free(p);
        return NULL;
    }
    if ((p[0][0][0]=(void**)malloc(n5*n4*n3*n2*sizeof(void*)))==NULL) {
        free(p[0][0]);
        free(p[0]);
        free(p);
        return NULL;
    }
    if ((p[0][0][0][0]=(void*)malloc(n5*n4*n3*n2*n1*size))==NULL) {
        free(p[0][0][0]);
        free(p[0][0]);
        free(p[0]);
        free(p);
        return NULL;
    }
    for (i5=0; i5<n5; i5++) {
        p[i5] = p[0]+n4*i5;
        for (i4=0; i4<n4; i4++) {
            p[i5][i4] = p[0][0]+n3*(i4+n4*i5);
            for (i3=0; i3<n3; i3++) {
                p[i5][i4][i3] = p[0][0][0]+n2*(i3+n3*(i4+n4*i5));
                for (i2=0; i2<n2; i2++)
		    p[i5][i4][i3][i2] = (char*)p[0][0][0][0]+
			size*n1*(i2+n2*(i3+n3*(i4+n4*i5)));
            }
        }
    }
    return p;
}

/* free a 5-d array */
void free5 (void *****p)
{
    free(p[0][0][0][0]);
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/* allocate a 6-d array */
void ******alloc6 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6,
		   size_t size)
{
    size_t i6,i5,i4,i3,i2;
    void ******p;

    if ((p=(void******)malloc(n6*sizeof(void*****)))==NULL)
        return NULL;

    if ((p[0]=(void*****)malloc(n6*n5*sizeof(void****)))==NULL) {
	free(p);
	return NULL;
    }

    if ((p[0][0]=(void****)malloc(n6*n5*n4*sizeof(void***)))==NULL) {
	free(p[0]);
	free(p);
	return NULL;
    }
    if ((p[0][0][0]=(void***)malloc(n6*n5*n4*n3*sizeof(void**)))==NULL) {
	free(p[0][0]);
	free(p[0]);
	free(p);
	return NULL;
    }
    if ((p[0][0][0][0]=(void**)malloc(n6*n5*n4*n3*n2*sizeof(void*)))==NULL) {
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
	return NULL;
    }
    if ((p[0][0][0][0][0]=(void*)malloc(n6*n5*n4*n3*n2*n1*size))==NULL) {
	free(p[0][0][0][0]);
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
	return NULL;
    }

    for (i6=0; i6<n6; i6++) {
	p[i6] = p[0]+n5*i6;
	for (i5=0; i5<n5; i5++) {
	    p[i6][i5] = p[0][0]+n4*(i5+n5*i6);
	    for (i4=0; i4<n4; i4++) {
		p[i6][i5][i4] = p[0][0][0]+n3*(i4+n4*(i5+n5*i6));
		for (i3=0; i3<n3; i3++) {
		    p[i6][i5][i4][i3] = p[0][0][0][0]
                        +n2*(i3+n3*(i4+n4*(i5+n5*i6)));
		    for (i2=0; i2<n2; i2++)
                        p[i6][i5][i4][i3][i2] =
			    (char*)p[0][0][0][0][0]+
			    size*n1*(i2+n2*(i3+n3*(i4+n4*(i5+n5*i6))));
		}
	    }
	}
    }
    return p;
}

/* free a 6-d array */
void free6 (void ******p)
{
    free(p[0][0][0][0][0]);
    free(p[0][0][0][0]);
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/* Zero-initialized allocation functions */

/* allocate a 1-d array and initialize to zero */
void *alloc1_zero (size_t n1, size_t size)
{
    void *p;

    if ((p=calloc(n1, size))==NULL)
	return NULL;
    return p;
}

/* allocate a 2-d array and initialize to zero */
void **alloc2_zero (size_t n1, size_t n2, size_t size)
{
    size_t i2;
    void **p;

    if ((p=(void**)calloc(n2, sizeof(void*)))==NULL)
	return NULL;
    if ((p[0]=(void*)calloc(n2*n1, size))==NULL) {
        free(p);
        return NULL;
    }
    for (i2=0; i2<n2; i2++)
	p[i2] = (char*)p[0]+size*n1*i2;
    return p;
}

/* allocate a 3-d array and initialize to zero */
void ***alloc3_zero (size_t n1, size_t n2, size_t n3, size_t size)
{
    size_t i3,i2;
    void ***p;

    if ((p=(void***)calloc(n3, sizeof(void**)))==NULL)
	return NULL;
    if ((p[0]=(void**)calloc(n3*n2, sizeof(void*)))==NULL) {
        free(p);
        return NULL;
    }
    if ((p[0][0]=(void*)calloc(n3*n2*n1, size))==NULL) {
        free(p[0]);
        free(p);
        return NULL;
    }

    for (i3=0; i3<n3; i3++) {
        p[i3] = p[0]+n2*i3;
        for (i2=0; i2<n2; i2++)
	    p[i3][i2] = (char*)p[0][0]+size*n1*(i2+n2*i3);
    }
    return p;
}

/* allocate a 4-d array and initialize to zero */
void ****alloc4_zero (size_t n1, size_t n2, size_t n3, size_t n4, size_t size)
{
    size_t i4,i3,i2;
    void ****p;

    if ((p=(void****)calloc(n4, sizeof(void***)))==NULL)
	return NULL;
    if ((p[0]=(void***)calloc(n4*n3, sizeof(void**)))==NULL) {
        free(p);
        return NULL;
    }
    if ((p[0][0]=(void**)calloc(n4*n3*n2, sizeof(void*)))==NULL) {
        free(p[0]);
        free(p);
        return NULL;
    }
    if ((p[0][0][0]=(void*)calloc(n4*n3*n2*n1, size))==NULL) {
        free(p[0][0]);
        free(p[0]);
        free(p);
        return NULL;
    }
    for (i4=0; i4<n4; i4++) {
        p[i4] = p[0]+i4*n3;
        for (i3=0; i3<n3; i3++) {
            p[i4][i3] = p[0][0]+n2*(i3+n3*i4);
            for (i2=0; i2<n2; i2++)
		p[i4][i3][i2] = (char*)p[0][0][0]+
		    size*n1*(i2+n2*(i3+n3*i4));
        }
    }
    return p;
}

/* allocate a 5-d array and initialize to zero */
void *****alloc5_zero (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size)
{
    size_t i5,i4,i3,i2;
    void *****p;

    if ((p=(void*****)calloc(n5, sizeof(void****)))==NULL)
	return NULL;
    if ((p[0]=(void****)calloc(n5*n4, sizeof(void***)))==NULL) {
        free(p);
        return NULL;
    }
    if ((p[0][0]=(void***)calloc(n5*n4*n3, sizeof(void**)))==NULL) {
        free(p[0]);
        free(p);
        return NULL;
    }
    if ((p[0][0][0]=(void**)calloc(n5*n4*n3*n2, sizeof(void*)))==NULL) {
        free(p[0][0]);
        free(p[0]);
        free(p);
        return NULL;
    }
    if ((p[0][0][0][0]=(void*)calloc(n5*n4*n3*n2*n1, size))==NULL) {
        free(p[0][0][0]);
        free(p[0][0]);
        free(p[0]);
        free(p);
        return NULL;
    }
    for (i5=0; i5<n5; i5++) {
        p[i5] = p[0]+n4*i5;
        for (i4=0; i4<n4; i4++) {
            p[i5][i4] = p[0][0]+n3*(i4+n4*i5);
            for (i3=0; i3<n3; i3++) {
                p[i5][i4][i3] = p[0][0][0]+n2*(i3+n3*(i4+n4*i5));
                for (i2=0; i2<n2; i2++)
		    p[i5][i4][i3][i2] = (char*)p[0][0][0][0]+
			size*n1*(i2+n2*(i3+n3*(i4+n4*i5)));
            }
        }
    }
    return p;
}

/* allocate a 6-d array and initialize to zero */
void ******alloc6_zero (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6,
		   size_t size)
{
    size_t i6,i5,i4,i3,i2;
    void ******p;

    if ((p=(void******)calloc(n6, sizeof(void*****)))==NULL)
        return NULL;

    if ((p[0]=(void*****)calloc(n6*n5, sizeof(void****)))==NULL) {
	free(p);
	return NULL;
    }

    if ((p[0][0]=(void****)calloc(n6*n5*n4, sizeof(void***)))==NULL) {
	free(p[0]);
	free(p);
	return NULL;
    }
    if ((p[0][0][0]=(void***)calloc(n6*n5*n4*n3, sizeof(void**)))==NULL) {
	free(p[0][0]);
	free(p[0]);
	free(p);
	return NULL;
    }
    if ((p[0][0][0][0]=(void**)calloc(n6*n5*n4*n3*n2, sizeof(void*)))==NULL) {
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
	return NULL;
    }
    if ((p[0][0][0][0][0]=(void*)calloc(n6*n5*n4*n3*n2*n1, size))==NULL) {
	free(p[0][0][0][0]);
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
	return NULL;
    }

    for (i6=0; i6<n6; i6++) {
	p[i6] = p[0]+n5*i6;
	for (i5=0; i5<n5; i5++) {
	    p[i6][i5] = p[0][0]+n4*(i5+n5*i6);
	    for (i4=0; i4<n4; i4++) {
		p[i6][i5][i4] = p[0][0][0]+n3*(i4+n4*(i5+n5*i6));
		for (i3=0; i3<n3; i3++) {
		    p[i6][i5][i4][i3] = p[0][0][0][0]
                        +n2*(i3+n3*(i4+n4*(i5+n5*i6)));
		    for (i2=0; i2<n2; i2++)
                        p[i6][i5][i4][i3][i2] =
			    (char*)p[0][0][0][0][0]+
			    size*n1*(i2+n2*(i3+n3*(i4+n4*(i5+n5*i6))));
		}
	    }
	}
    }
    return p;
}

