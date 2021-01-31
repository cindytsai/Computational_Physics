/*-------------------------------------------------------------------------

    String.c                                  November, 2017

    OpenGL Animation

--------------------------------------------------------------------------*/


#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef MACOSX
#include <OpenGl/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

GLint    N     = 1000;      // the number of grid points in x-direction
GLdouble L     = 8.0;       // length of the string
GLdouble T     = 100.0;       // string tension
GLdouble mu    = 0.1;       // mass density
GLdouble g     = 9.8;       // gravity constant
GLdouble k     = 5000.0;    // String hook constant
GLdouble A0    = 1.0;       // amplitude of the wave

GLint    nl    = 4;         // number of half-wavelength (for Sinusoidal)
GLdouble X0    = 4.0;       // peak location (for Plucked and Gaussian)
GLdouble sigma = 0.5;       // half-width of the peak (for Gaussian)

// The following are computed in main.
GLdouble bound = 0.0;       // bounds of string amplitude
GLdouble c     = 0.0;       // wave velocity
GLdouble dx    = 0.0;       // space increment
GLdouble dt    = 0.0;       // time increment (for exact finite difference)
GLdouble epsilon = 0.0;

// Arrays for the string.
double *x;          // coordinate of x positions.
double *U0;         // string amplitudes
double *dU0;        // time-derivation evaluated at t=0

int countPrint = 0;

/*--------------------------------------------------------------------------
 *  Utility subroutine
 */
void *memalloc(size_t n){
    void *v;

    if ((v = malloc(n)) == NULL) {
        printf("!!! Not enough memory.\n");
        exit(1);
    }
    return v;
}

/*--------------------------------------------------------------------------
 *  Prepare the initial waves.
 */
void Sinusoidal(int N, double *x, double *U0, double *dU0, double A0, double L, double dx, int nl){
    int  i;

    for (i=0; i <= N; i++) {
        x[i]   = i * dx;
        U0[i]  = A0*sin(M_PI*nl/L * i * dx);
        dU0[i] = 0.0;
    }
}

void Plucked(int N, double *x, double *U0, double *dU0, double A0, double L, double dx, double X0){
    int  i;

    for (i=0; i <= N; i++) {
        x[i]   = i * dx;
        U0[i]  = (x[i] < X0) ? x[i]*A0/X0 : (L-x[i])*A0/(L-X0);
        dU0[i] = 0.0;
    }
}

void Gaussian(int N, double *x, double *U0, double *dU0, double A0, double L, double dx, double X0, double sigma){
    double a;
    int    i;

    for (i=0; i <= N; i++) {
        x[i]   = i * dx;
        a      = (x[i]-X0)/sigma;
        U0[i]  = A0*exp(-0.5*a*a);
        dU0[i] = 0.0;
    }
}

/*--------------------------------------------------------------------------
 * Solve PDE of vibrating string, with gravity
 */
void PDE(int cnt, int N, double *U0, double *dU0, double epsilon, double dt, double bound){
    static double *U, *Uold;
    int    j;

    if (U    == NULL) U    = memalloc((N+1)*sizeof(double));
    if (Uold == NULL) Uold = memalloc((N+1)*sizeof(double));

    U[0] = 0.0;
    U[N] = 0.0;

    if (cnt == 0) {
        for (j=0; j <= N; j++){
            U[j] = 0.5 * epsilon * (U0[j+1] + U0[j-1]) + (1.0 - epsilon) * U0[j] + dt * dU0[j] - 0.5 * g * pow(dt, 2);
        }
    }
    else {
        for (j=0; j <= N; j++) {
            U[j] = epsilon * (U0[j+1] + U0[j-1]) + 2.0 * (1.0 - epsilon) * U0[j] - Uold[j] - g * pow(dt, 2);
            // if (fabs(U[j]) > bound) {
            //     printf("String out of range: %d: %f\n", j, U[j]);
            //     exit(1);
            // }
        }
    }

    // Uold   -> u_j-1
    // U0     -> u_j
    // U      -> u_j+1 , what we are calculating
    memcpy(Uold, U0, sizeof(double)*(N + 1));
    memcpy(U0,   U,  sizeof(double)*(N + 1));

    // Print out maximum in U0
    // if(countPrint <= 10000){
    //     GLdouble maxU0 = -1000.0;
    //     FILE *output;
    //     output = fopen("g_gaussian_maxA0.txt", "a");
    //     for(int l = 0; l <= N; l = l+1){
    //         if(U0[l] > maxU0){
    //             maxU0 = U0[l];
    //         }
    //     }
    //     fprintf(output, "%.5lf\n", maxU0);
    //     countPrint = countPrint + 1;
    //     fclose(output);
    // }
    // else{
    //     exit(0);
    // }
}

/*--------------------------------------------------------------------------
 * Handler for window-repaint event. Call back when the window first appears
 * and whenever the window needs to be re-painted.
 */
void display(void){
    static int cnt;
    double x1, x2;
    int i;

    glClearColor(1.0, 1.0, 1.0, 1.0);                    // white background
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  // clear the color buffer
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,bound,0,0,0,0,1,0);      // 3d view direction
    
    x1 = -L*0.10;
    x2 =  L*1.10;
    
    glBegin(GL_LINES);          // Draw x-axis
    glColor3f(0.7,0.7,0.7);
    glVertex2f(x1,0);
    glVertex2f(x2,0);
    glEnd();
    
    glBegin(GL_LINES);          // Draw y-axis
    glColor3f(0.7,0.7,0.7);
    glVertex2f(0,-bound);
    glVertex2f(0, bound);
    glEnd();

    glBegin(GL_LINES);
    glColor3f(0.7,0.7,0.7);
    glVertex2f(x[N],-bound);
    glVertex2f(x[N], bound);
    glEnd();

    glBegin(GL_LINES);
    glColor3f(0.0,0.6,0.0);
    for(i=0; i < N; i++) {
        glVertex2f(x[i],   U0[i]);
        glVertex2f(x[i+1], U0[i+1]);
    }
    glEnd();
 
    glutSwapBuffers();    // Double buffered - swap the front and back buffers

    PDE(cnt, N, U0, dU0, epsilon, dt, bound);
    cnt = cnt + 1;
}

/*----------------------------------------------------------------------------
 * Handler for window re-size event. Called back when the window first appears
 * and whenever the window is re-sized with its new width and height
 */
void WindowSize(int w, int h){
    double x1, x2;

    x1 = -0.10*L;
    x2 =  1.10*L;

    glViewport(0, 0, w, h);       // initial draw position and the size
    glMatrixMode(GL_PROJECTION);    // project 3d into 2d
    glLoadIdentity();       // clean up the matrix data with identity
    glOrtho(x1,x2,-bound,bound,-bound,bound);
                // matrix box of (x1,x2),(y1,y2),(z1,z2)
    glMatrixMode(GL_MODELVIEW);     // matrix for rotation, scaling and translation
    glLoadIdentity();
}

/* Called back when there is no other event to be handled */
void idle(void){
    glutPostRedisplay();        // post a re-paint request to activate display()
}


int main(int argc, char** argv){
    int type = 1;

    glutInit(&argc, argv);                         // initialize GLUT

    printf("Please select the type of string:\n");
    printf(" 1. Sinusoidal wave.\n");
    printf(" 2. Plucked string.\n");
    printf(" 3. Gaussian wave.\n");
    printf("Your selection: ");
    scanf("%d", &type);
    fflush(stdout);

    x   = memalloc((N + 1) * sizeof(double));
    U0  = memalloc((N + 1) * sizeof(double));
    dU0 = memalloc((N + 1) * sizeof(double));
    c   = sqrt(T/mu);
    dx  = L / N;
    dt  = dx / c;
    bound   = A0 * 1.25;                             // Bounds of string amplitude
    epsilon = pow((dt * c / dx), 2);

    if (type == 1){
        Sinusoidal(N, x, U0, dU0, A0, L, dx, nl);
    }
    
    else if (type == 2){
        Plucked(N, x, U0, dU0, A0, L, dx, X0);
    }
    else{
        Gaussian(N, x, U0, dU0, A0, L, dx, X0, sigma);
    }

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);   // enable double buffered mode
    glutInitWindowSize(500,250);                     // window size
    glutInitWindowPosition(300,100);             // window position (top-left)
    glutCreateWindow("Vibrating String");          // creat window with title
    glutReshapeFunc(WindowSize);                   // callback for re-size event
    glutDisplayFunc(display);                      // callback for re-paint event
    glutIdleFunc(idle);                            // callback if no other event
    glutMainLoop();                                // enter event-processing loop
 
    return 0;
} 
