#include <stdio.h> 
#include <stdlib.h> 
#include <GL/glut.h>
#include <math.h>
#include <string.h>

// Basic Settings
// N    -> Length of the box
// R    -> Radius of the plate
// U    -> array to store the new calculate temperature
// U0   -> array of temperature we display
// dt   -> basic unit of time
// dx   -> basic unit of cordinate
// c    -> Propergating speed
int N = 100;
double R = 40;
double *U0;
double *U;
double dt, dx, c;

// For mapping temperature to colormap Rainbow
GLdouble maxTemp = 380.0;
GLdouble minTemp = 0.0;

struct RGBColor{
    GLfloat red;
    GLfloat green;
    GLfloat blue;
};

struct RGBColor tempTocolormap(double temperature){
	// Make temperature into a form 
	// like (0 ~ 1) -> (min ~ man) -> (Blue ~ Red)
	double temp01 = (temperature - minTemp) / (maxTemp - minTemp);
	struct RGBColor colormap = {.red = 0.0, .green = 0.0, .blue = 0.0};

	// Invert it due to the alogrithm
	temp01 = 1.0 - temp01;

	if(temp01 >= 0.0 && temp01 < 0.25){
		colormap.red = 255.0;
		colormap.green = 255.0 * (temp01 / 0.25);
		colormap.blue = 0.0;
	}
	if(temp01 >= 0.25 && temp01 < 0.5){
		colormap.red = 255.0 * (1 - ((temp01 - 0.25) / 0.25));
		colormap.green = 255.0;
		colormap.blue = 0.0;
	}
	if(temp01 >= 0.5 && temp01 < 0.75){
		colormap.red = 0.0;
		colormap.green = 255.0;
		colormap.blue = 255.0 * (temp01 - 0.5) / 0.25;
	}
	if(temp01 >= 0.75 && temp01 <= 1.0){
		colormap.red = 0.0;
		colormap.green = 255.0 * (1 - ((temp01 - 0.75) / 0.25));
		colormap.blue = 255.0;
	}

	return colormap;
}

void *memalloc(size_t n){
	void *v;
	if ((v = malloc(n)) == NULL) {
		printf("!!! Not enough memory.\n");
		exit(1);
	}
	return v;
}

double boundaryValue(int x, int y, int N){
	if((x - N/2) >= 0 && (y - N/2) >= 0){
		return 373.0;
	}
	else{
		return 273.0;
	}
}

int isOutsideBoundary(double x, double y, int N){
	if((pow(x - N/2, 2) + pow(y - N/2, 2)) > pow(R, 2)){
		return 1;
	} 
	else{
		return 0;
	}
}

int isInsideBoundary(double x, double y, int N){
	if(pow(x - N/2, 2) + pow(y - N/2, 2) < pow(R, 2)){
		return 1;
	}
	else{
		return 0;
	}
}

void heatDiffusion(){
	// Settings
	// fw, bw     -> Forward and Backward of the index
	// eps1~4     -> Epsilon for boundary {right, left, top, down}
	// f1~4       -> Value for {right, left, top, down}
	// corx,cory  -> The coordinate we are calculate
	// neix, neiy -> The neighborhood's coordinate
	// temp       -> Store temperary stuff when calculating epsilon
	int fw[N], bw[N];
	double eps1, eps2, eps3, eps4;
	double f1, f2, f3, f4;
	double corx, cory;
	double neix, neiy;
	double temp, sum;

	for(int i = 0; i < N; i = i+1) {
    	fw[i] = (i + 1) % N;
    	bw[i] = (i - 1 + N) % N;
  	}

  	// Calculate the new temperature array U[N*N]
  	for(int i = 0; i < N; i = i+1){
  		for(int j = 0; j < N; j = j+1){

  			// Initialize
  			temp = 0.0;
  			sum = 0.0;
  			
  			// Set the coordinate corx, cory
  			corx = (double)i + 0.5;
  			cory = (double)j + 0.5;

  			if(isInsideBoundary(corx, cory, N)){
  				// Check the right point
  				neix = (double)fw[i] + 0.5;
  				neiy = (double)j + 0.5;

  				if(isOutsideBoundary(neix, neiy, N)){
  					f1 = boundaryValue(fw[i], j, N);
  					temp = (double)N/2.0 + sqrt(pow(R, 2) - pow(cory - (double)N/2.0, 2));
  					eps1 = (temp - corx) * dx;
  				}
  				else{
  					f1 = U0[fw[i]+j*N];
  					eps1 = dx;
  				}

  				// Check the left point
  				neix = (double)bw[i] + 0.5;
  				neiy = (double)j + 0.5;

  				if(isOutsideBoundary(neix, neiy, N)){
  					f2 = boundaryValue(bw[i], j, N);
  					temp = (double)N/2.0 - sqrt(pow(R, 2) - pow(cory - (double)N/2.0, 2));
  					eps2 = (corx - temp) * dx;
  				}
  				else{
  					f2 = U0[bw[i]+j*N];
  					eps2 = dx;
  				}

  				// Check the top point
  				neix = (double)i + 0.5;
  				neiy = (double)fw[j] + 0.5;

  				if(isOutsideBoundary(neix, neiy, N)){
  					f3 = boundaryValue(i, fw[j], N);
  					temp = (double)N/2.0 + sqrt(pow(R, 2) - pow(corx - (double)N/2.0, 2));
  					eps3 = (temp - cory) * dx;
  				}
  				else{
  					f3 = U0[i+N*fw[j]];
  					eps3 = dx;
  				}

  				// Check the bottom point
  				neix = (double)i + 0.5;
  				neiy = (double)bw[j] + 0.5;

  				if(isOutsideBoundary(neix, neiy, N)){
  					f4 = boundaryValue(i, bw[j], N);
  					temp = (double)N/2.0 - sqrt(pow(R, 2) - pow(corx - (double)N/2.0, 2));
  					eps4 = (cory - temp) * dx;
  				}
  				else{
  					f4 = U0[i+N*bw[j]];
  					eps4 = dx;
  				}

  				// Update the temperature inside the boundary
  				sum = (2.0 / (eps1 * eps2 * (eps1 + eps2))) * (eps2 * f1 + eps1 * f2 - (eps1 + eps2) * U0[i+j*N]);
  				sum = sum + (2.0 / (eps3 * eps4 * (eps3 + eps4))) * (eps3 * f4 + eps4 * f3 - (eps3 + eps4) * U0[i+j*N]);
  				U[i+j*N] = pow(c, 2) * dt * sum + U0[i+j*N];
  			}

  		}
  	}

  	// Copy memery in U to U0, since we display U0
  	memcpy(U0, U, sizeof(double) * (N*N));
}


/*--------------------------------------------------------------------------
 * Handler for window-repaint event. Call back when the window first appears
 * and whenever the window needs to be re-painted.
 */
void display(void){
	// Settings
	// cnt         -> For time series
	// color_temp  -> Hold the colormap red, green blue
	static int cnt;
	struct RGBColor color_temp;
  	
  	glClearColor(1.0, 1.0, 1.0, 1.0);                    // white background
  	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  // clear the color buffer
	
  	glMatrixMode(GL_MODELVIEW);
  	glLoadIdentity();
  	gluLookAt(0, 0, 1.0, 0, 0, 0, 0, 1, 0);      		 // 3d view direction
	
	// Start from here, draw a box with length N
  	glBegin(GL_LINES);			
	glColor3f(0.7,0.7,0.7);
	glVertex2f(0, 0);
	glVertex2f(N, 0);
  	glEnd();
  	
  	glBegin(GL_LINES);
	glColor3f(0.7,0.7,0.7);
	glVertex2f(N, 0);
	glVertex2f(N, N);
  	glEnd();

  	glBegin(GL_LINES);
  	glColor3f(0.7,0.7,0.7);
  	glVertex2f(N, N);
  	glVertex2f(0, N);
  	glEnd();

  	glBegin(GL_LINES);
  	glColor3f(0.7,0.7,0.7);
  	glVertex2f(0, N);
  	glVertex2f(0, 0);
  	glEnd();

  	// Plot the temperature colormap
  	for(int i = 0; i < N; i = i+1){
  		for(int j = 0; j < N; j = j+1){
  			// Get color
  			color_temp = tempTocolormap(U0[i+j*N]);

  			// Draw
  			glBegin(GL_QUADS);
  			glColor3ub(color_temp.red, color_temp.green, color_temp.blue);
  			glNormal3f(0.0, 0.0, 1.0);
  			glVertex3f(  i,   j, 0);
  			glVertex3f(i+1,   j, 0);
  			glVertex3f(i+1, j+1, 0);
  			glVertex3f(  i, j+1, 0);
  			glEnd();
  		}
  	}

  	// Plot the fixed temperature ring, 
  	// the upper right temperature      -> 373.0
  	// the rest of the ring temperature -> 273.0
  	struct RGBColor color373, color273;
  	color373 = tempTocolormap(373.0);
  	color273 = tempTocolormap(273.0);

  	for(double angle = 0.0; angle <= 2*M_PI; angle = angle + 0.01){
  		glBegin(GL_LINES);
  		if(angle >= 0.0 && angle < M_PI / 2.0){
  			glColor3ub(color373.red, color373.green, color373.blue);
  		}
  		else{
  			glColor3ub(color273.red, color273.green, color273.blue);
  		}
  		glLineWidth(10.0);
  		glVertex2f(R * cos(angle) + (double)N/2.0, R * sin(angle) + (double)N/2.0);
  		glVertex2f(R * cos(angle+0.01) + (double)N/2.0, R * sin(angle+0.01) + (double)N/2.0);
  		glEnd();
  	}
  	
  	glutSwapBuffers();	  // Double buffered - swap the front and back buffers

  	heatDiffusion();

  	cnt = cnt + 1;
}

/*----------------------------------------------------------------------------
 * Handler for window re-size event. Called back when the window first appears
 * and whenever the window is re-sized with its new width and height
 */
void WindowSize(int w, int h){
  	double x1, x2;

  	x1 = - 0.05 * (double)N;
  	x2 = 1.05 * (double)N;

  	glViewport(0, 0, w, h);      				// initial draw position and the size
  	glMatrixMode(GL_PROJECTION);  				// project 3d into 2d
  	glLoadIdentity();			  				// clean up the matrix data with identity
  	glOrtho(x1, x2, x1, x2, -1.0, 1.0);			// matrix box of (x1,x2),(y1,y2),(z1,z2)
  	glMatrixMode(GL_MODELVIEW); 				// matrix for rotation, scaling and translation
  	glLoadIdentity();
}

/* Called back when there is no other event to be handled */
void idle(void){
  	glutPostRedisplay();        // post a re-paint request to activate display()
}

int main(int argc, char** argv){
	
	// Initialize
	U0 = memalloc((N * N) * sizeof(double));
	U = memalloc((N * N) * sizeof(double));

	dt = 1.0 / (double) N;
	dx = 1.0;
	c = 1.0;

	// Assign the plate to specific environment
	for(int i = 0; i < N; i = i+1){
		for(int j = 0; j < N; j = j+1){
			if(isOutsideBoundary((double)i+0.5, (double)j+0.5, N)){
				U0[i+j*N] = -1.0;
			}
			else{
				U0[i+j*N] = 0.0;
			}
		}
	}
	
	// GLUT
	glutInit(&argc, argv);                         // initialize GLUT
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);   // enable double buffered mode
  	glutInitWindowSize(1000, 1000);         	   // window size
  	glutInitWindowPosition(300,100);         	   // window position (top-left)
  	glutCreateWindow("Heat Diffusion Simulation"); // creat window with title
  	glutReshapeFunc(WindowSize);                   // callback for re-size event
  	glutDisplayFunc(display);                      // callback for re-paint event
  	glutIdleFunc(idle);                            // callback if no other event
  	glutMainLoop();                                // enter event-processing loop


	return 0;
}