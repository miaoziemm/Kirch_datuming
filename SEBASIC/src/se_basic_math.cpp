#include "../include/se_basic_math.h"

/** Gaussian ("normally") distributed random values with mean 0.0 and
    standard deviation 1.0.

    This uses the polar method of G. E. P. Box, M. E. Muller, and G. Marsaglia,
    as described by Donald E. Knuth in The Art of Computer Programming, 
    Volume 2: Seminumerical Algoztrhms, section 3.4.1, subsection C,
    algoztrhm P. 
*/
real rand_normal()
{
    static unsigned char have_next = 0;
    static real next;
    real v1, v2, s, tmp;

    if ( have_next ) {
        have_next = 0;
        return next;
    }
    
    do { 
        v1 = 2 * rand1() - 1;
        v2 = 2 * rand1() - 1;
        s = v1*v1 + v2*v2;
    } while ( s >= 1 );

    if( ! FISZERO(s) ) {
        tmp = sqrt( -2 * log(s) / s );
        next = v2 * tmp;
        have_next = 1;
        return v1 * tmp;
    }
    else {
        return 0;
    }
}


reala percentile(real p, reala* samples, size_t ns)
{
    ssize_t k, lo, hi;

    k = (ssize_t) (ns * p / 100.0);
    if(k < 0)   k = 0;
    if(k >= (ssize_t)ns) k = ns - 1;

    for (lo = 0, hi = ns - 1; lo < hi; ) {
        reala xk = samples[k];
        ssize_t i = lo;
        ssize_t j = hi;
        do {
            while (samples[i] < xk) {
                i++;
            }
            while (samples[j] > xk) {
                j--;
            }
            if (i <= j) {
                reala tmp = samples[i];
                samples[i] = samples[j];
                samples[j] = tmp;
                i++;
                j--;
            } else { // spare a comparison
                break;
            }
        } while (i <= j);

        if (j < k) {
            lo = i;
        }
        if (k < i) {
            hi = j;
        }
    }

    return samples[k];
}
// vectors -> vector ops
void v3d_cross(v3d_t *u, const v3d_t *v, const v3d_t *w)
{
        
    u->v[0] =   v->v[1] * w->v[2] - v->v[2] * w->v[1];
    u->v[1] = - v->v[0] * w->v[2] + v->v[2] * w->v[0];
    u->v[2] =   v->v[0] * w->v[1] - v->v[1] * w->v[0];

    if (FISZERO(v3d_norm(u))) {

        if ( v3d_norm(v) > v3d_norm(w) ) {
            v3d_get_orthogonal(u,v);
        } else {
            v3d_get_orthogonal(u,w);            
        }
    }
}

void v3d_get_orthogonal(v3d_t *u, const v3d_t *v)
{
    int min_ind=1;
    double min=fabs(v->v[0]);
    if (min > fabs(v->v[1])) {
        min = fabs(v->v[1]);
        min_ind=2;
    }
    if (min > fabs(v->v[2])) {
        min = fabs(v->v[2]);
        min_ind=3;
    }
    
    switch(min_ind) {
    case 1:
        u->v[0]=0;
        u->v[1]= v->v[2];
        u->v[2]=-v->v[1];
        break;
    case 2:
        u->v[0]= v->v[2];
        u->v[1]=0;
        u->v[2]=-v->v[0];
        break;
    case 3:
        u->v[0]= v->v[1];
        u->v[1]=-v->v[0];
        u->v[2]=0;
        break;
    default:
        u->v[0]=0;
        u->v[1]=0;
        u->v[2]=0;
        break;
    } 
}

void v3d_get_orthogonal_pair(v3d_t *u, v3d_t *w, const v3d_t *v)
{
    v3d_get_orthogonal(u, v);
    v3d_cross(w, v, u);
}

// Matrix ops
int m3d_is_sym(const m3d_t *A)
{

    double diff;

    diff = ( FISZERO(fabs(A->m[1])) && FISZERO(fabs(A->m[3])) ) ? 0 : 
        ( A->m[1] - A->m[3] ) / ( fabs(A->m[1]) + fabs(A->m[3]) ) ; 
    if (fabs(diff)  > 1e-6) {
        INFOV((10, "Matrix is not symmetric entry (1,2) != (2,1)"));
        return 0;
    }

    diff = ( FISZERO(fabs(A->m[2])) && FISZERO(fabs(A->m[6])) ) ? 0 : 
        ( A->m[2] - A->m[6] ) / ( fabs(A->m[2]) + fabs(A->m[6]) ) ; 
    if (fabs(diff)  > 1e-6) {
        INFOV((10, "Matrix is not symmetric entry (1,3) != (3,1)"));
        return 0;
    }

    diff = ( FISZERO(fabs(A->m[5])) && FISZERO(fabs(A->m[7])) ) ? 0 : 
        ( A->m[5] - A->m[7] ) / ( fabs(A->m[5]) + fabs(A->m[7]) ) ; 
    if (fabs(diff)  > 1e-6) {
        INFOV((10, "Matrix is not symmetric entry (2,3) != (3,2)"));
        return 0;
    }

    return 1;

}



// Matrix -> Matrix ops

static double get_2d_det(double a11, double a12, 
                              double a21, double a22) {
    return a11*a22 - a12*a21;
}

int m3d_inverse(m3d_t *A_inv, const m3d_t *A)
{

    double d=m3d_determinant(A);

    if (FISZERO(d)) {

        m3d_init_zero(A_inv);
        return 1;

    } else {

        d = 1/d;

        A_inv->m[0] = get_2d_det(A->m[4],A->m[5],A->m[7],A->m[8]) * d;
        A_inv->m[1] = get_2d_det(A->m[2],A->m[1],A->m[8],A->m[7]) * d;
        A_inv->m[2] = get_2d_det(A->m[1],A->m[2],A->m[4],A->m[5]) * d;

        A_inv->m[3] = get_2d_det(A->m[5],A->m[3],A->m[8],A->m[6]) * d;
        A_inv->m[4] = get_2d_det(A->m[0],A->m[2],A->m[6],A->m[8]) * d;
        A_inv->m[5] = get_2d_det(A->m[2],A->m[0],A->m[5],A->m[3]) * d;

        A_inv->m[6] = get_2d_det(A->m[3],A->m[4],A->m[6],A->m[7]) * d;
        A_inv->m[7] = get_2d_det(A->m[1],A->m[0],A->m[7],A->m[6]) * d;
        A_inv->m[8] = get_2d_det(A->m[0],A->m[1],A->m[3],A->m[4]) * d;

        return 0;
    }
}


/**
 * for the 3x3 symmetric case use the characteristic polynomial 
 * det(m*I - M) = m^3 - m^2 tr(M) - .5m (tr(M^2)-tr(M)^2) - det(M)
 * use the tranform M=f B + g I, M and B have the same eigenvectors  
 * and eigenvalues of M are eigenvalues of B times p plus q
 * use g = tr(M)/3 and f = tr((A-g I)^2/6)^.5
 * then 
 * det (b*I+B) = b^3 - 3b - det(B)
 * solve using b=2cos(3t)
 * get b = 2cos(1/3 (arccos(det(B)/2) + 2 pi n)) for n=0,1,2
 *
 */
void m3d_sym_eigenvals(const m3d_t *M, double *ev) 
{
    double f, g, h, r, t;
    m3d_t B;
    
    h = M->m[1]*M->m[1] + M->m[2]*M->m[2] + M->m[5]*M->m[5];
    g = m3d_trace(M)/3;
 
    // check if M is diagonal
    if ( (h/fabs(g))<1e-6 ) {
        ev[0] = M->m[0];
        ev[1] = M->m[4];
        ev[2] = M->m[8];
        return;
    }

    f = sqrt( ( (M->m[0]-g)*(M->m[0]-g) + (M->m[4]-g)*(M->m[4]-g) + (M->m[8]-g)*(M->m[8]-g) + 2*h ) / 6 );
    
    B = *M;
    B.m[0] -= g;
    B.m[4] -= g;
    B.m[8] -= g;
    
    r = m3d_determinant(&B) / (2*f*f*f);

    if (r<=-1) {
        t = M_PI/3.0;
    } else if (r>=1) {
        t = 0;
    } else {
        t = acos(r)/3.0;
    }
    
    ev[0] = g + 2*f*cos(t);
    ev[1] = g + 2*f*cos(t + 2.0*M_PI/3.0);
    ev[2] = 3*g - ev[0] - ev[1];

}

static double m4d_sym_find_largest_off_entry(const m4d_t *M, int *m_i, int *m_j)
{
    int i,j;
    double ent_norm;

    (*m_i)=0;
    (*m_j)=1;
    ent_norm = fabs(M->m[(*m_i)*4+(*m_j)]);

    for (i=0;i<4;i++) {
        for (j=i+1; j<4; j++) {
            if (ent_norm < fabs(M->m[i*4+j]) ) {
                ent_norm = fabs(M->m[i*4+j]);
                (*m_i) = i;
                (*m_j) = j;
            }
        }
    }
    
    return ent_norm;
}

static double m4d_sym_find_largest_diag_entry(const m4d_t *M) 
{
    int i;
    double ent_norm;

    ent_norm=0;
    for (i=0;i<4;i++) {
        if (ent_norm < fabs(M->m[4*i+i])) {
            ent_norm = fabs(M->m[4*i+i]);
        }
    }

    return ent_norm;
}

/**
 * Use Jacobi eigenvalue algorithm: Apply Givens rotations to zero
 * out the non-diagonal entries -- not possible to completely zero
 * out the off-diagonal entries, but this algorithm boosts the size
 * of the diagonal entries compared to the off-diagonal entries.
 * The diagonal entries will be close to the eigenvalues, after several
 * iterations. The same rotations will convert the identity matrix
 * to a matrix that contains the eigenvectors as columns.
 * 
 */
void m4d_sym_eigensystem(const m4d_t *M, double *ev, m4d_t *E)
{
    int i,j,iter;
    m4d_t S, Sp;

    if (E) m4d_diag_init(E, 1, 1, 1, 1);

    memcpy(S.m, M->m, 16*sizeof(double));
    iter=0;
    while ( ( m4d_sym_find_largest_off_entry(&S, &i, &j) 
              > 1e-8 * m4d_sym_find_largest_diag_entry(&S)) 
            && ( iter<100 ) ) {

        int k;
        double th, s, c;

        memcpy(Sp.m, S.m, 16*sizeof(double));

        if (FISZERO(S.m[4*i+i]-S.m[4*j+j])) {
            th = M_PI/4;
        } else {
            th = .5*atan(2*S.m[4*i+j]/(S.m[4*j+j]-S.m[4*i+i]));
        }

        sincos(th,&s,&c);

        Sp.m[4*i+i]=c*c*S.m[4*i+i] - 2*s*c*S.m[4*i+j] + s*s*S.m[4*j+j];
        Sp.m[4*j+j]=s*s*S.m[4*i+i] + 2*s*c*S.m[4*i+j] + c*c*S.m[4*j+j];

        Sp.m[4*i+j]=(c*c-s*s)*S.m[4*i+j] + s*c*(S.m[4*i+i] - S.m[4*j+j]);
        Sp.m[4*j+i] = Sp.m[4*i+j];
        
        for (k=0; k<4; k++) {
           
            if (E) {
                double Eik, Ejk;

                Eik = c*E->m[4*k+i] - s*E->m[4*k+j];
                Ejk = s*E->m[4*k+i] + c*E->m[4*k+j];

                E->m[4*k+i] = Eik;
                E->m[4*k+j] = Ejk;
            }

            if ((k==i)||(k==j)) continue;

            Sp.m[4*i+k] = c*S.m[4*i+k] - s*S.m[4*j+k];
            Sp.m[4*k+i] = Sp.m[4*i+k];

            Sp.m[4*j+k] = s*S.m[4*i+k] + c*S.m[4*j+k];
            Sp.m[4*k+j] = Sp.m[4*j+k];
                
        }

        memcpy(S.m, Sp.m, 16*sizeof(double));
        iter++;

    }

    ev[0] = S.m[ 0];
    ev[1] = S.m[ 5];
    ev[2] = S.m[10];
    ev[3] = S.m[15];


}
