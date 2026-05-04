#include "../include/se_geometry.h"


int se_find_closest_pt_on_planes(v3d_t *b, 
                                   const v3d_t *p1, const v3d_t *n1,
                                   const v3d_t *p2, const v3d_t *n2)
{
    v3d_t nv; 
    double n1n1, n2n2, nvnv, n1n2; 
    double p1n1, p2n2;
    double c1, c2, det;

    v3d_cross(&nv, n1, n2);  //叉乘
    nvnv = v3d_dot(&nv,&nv);  //点乘 nv*nv

    if (FISZERO(nvnv)) return 1;

    n1n1 = v3d_dot(n1,n1);
    n2n2 = v3d_dot(n2,n2);

    n1n2 = v3d_dot(n1,n2);
    
    det = n1n1*n2n2 - n1n2*n1n2;
    
    if (FISZERO(det)) return 1;

    p1n1 = v3d_dot(p1,n1);
    p2n2 = v3d_dot(p2,n2);

    c1 = (n2n2*p1n1 - n1n2*p2n2)/det;
    c2 = (n1n1*p2n2 - n1n2*p1n1)/det;

    v3d_assign_scaled(b, c1, n1);
    v3d_accum_scaled(b, c2, n2);

    v3d_accum_scaled(b, 
                             ( v3d_dot(p1,&nv) + v3d_dot(p2,&nv) ) / ( 2*nvnv ),
                             &nv);

    return 0;
 
}


static double get_step(const v3d_t * grad, const v3d_t *n, const m3d_t *M, const v3d_t *b_m_p)
{

    double A, B, C;

    A = .5 * v3d_m3d_v3d_inner(grad, M, grad) ;

    B = v3d_dot(n, grad) + v3d_m3d_v3d_inner(grad, M, b_m_p) ;    
    
    C = v3d_dot(n, b_m_p) + .5*v3d_m3d_sym_inner(M, b_m_p) ; 

    
    if (FISZERO(A)) {
        if (FISZERO(B)) { 
            return 0;
        } else {
            return -C/B;
        }
    } else {
        if ((B*B-4*A*C)<0) {
            return -B/(2*A);
        } else {
            return (-B + sqrt(B*B-4*A*C))/(2*A);
        }
    }
        
}

int se_find_closest_pt_on_paraboloids(v3d_t *b, double *F1, double *F2, 
                                        const v3d_t *p1, const v3d_t *n1, const m3d_t *M1,
                                        const v3d_t *p2, const v3d_t *n2, const m3d_t *M2,
                                        double tol)
{

    const int max_iter_count=10;

    if (se_find_closest_pt_on_planes(b, p1, n1, p2, n2)) {
        v3d_add_scaled(b, .5, p1, .5, p2);
    }

    int done=0;
    int it=0;
    do {
        v3d_t b_m_p1, b_m_p2, tmp, g1, g2, b1,b2;

        v3d_subtract(&b_m_p1, b, p1);  //b_m_p1=b-p1;
        v3d_subtract(&b_m_p2, b, p2);  //b_m_p2=b-p2;

        *F1 = se_eval_quadratic_form(NULL, 0, n1, M1, &b_m_p1);
        *F2 = se_eval_quadratic_form(NULL, 0, n2, M2, &b_m_p2);

        if ( ( (fabs(*F1)+fabs(*F2)) < tol) || (it>max_iter_count) ) {
            done = 1;
        } else {
            m3d_v3d_mult(&tmp, M1, &b_m_p1);
            v3d_add(&g1, n1, &tmp);

            m3d_v3d_mult(&tmp, M2, &b_m_p2);
            v3d_add(&g2, n2, &tmp);

            v3d_add_scaled(&b1,1,b,get_step(&g1, n1, M1, &b_m_p1),&g1);

            v3d_add_scaled(&b2,1,b,get_step(&g2, n2, M2, &b_m_p2),&g2);

            v3d_add_scaled(b,.5,&b1,.5,&b2);
        }
        it++;

    } while (!done);


    return 0;


}
