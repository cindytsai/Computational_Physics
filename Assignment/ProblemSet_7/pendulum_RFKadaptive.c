#include <stdio.h> 
#include <stdlib.h> 
#include <GL/glut.h>
#include <math.h>

GLdouble mass = 0.5;         // particle mass
GLdouble g  = 9.8;           // gravitational constant
GLdouble L  = 6.0;           // string length

GLdouble E;                  // total energy of the system
GLdouble E0;                 // Initial total energy of the system
GLdouble theta;              // initial angle
GLdouble W;                  // initial angular velocity

GLdouble delta0 = 1e-7;
GLdouble dt0 = 0.1;


/* Fourth Order Runge-Kutta Algorithm */
void RKF(GLdouble dt){
	double theta0 = theta;
	double W0 = W;
	double k1[4];
	double k2[4];
		
	for(int i = 0; i < 3; i = i+1){
		k1[i] = dt * (-(g/L) * sin(theta));
		k2[i] = dt * W;
		W  = W0 + 0.5 * k1[i];
		theta = theta0 + 0.5 * k2[i];
	}
	W = W0 + k1[2];
	theta = theta0 + k2[2];
	k1[3] = dt * (-(g/L) * sin(theta));
	k2[3] = dt * W;
	
	// Calulate W, and theta
	W = W0 + (k1[0] + 2.0 * k1[1] + 2.0 * k1[2] + k1[3]) / 6.0;
	theta = theta0 + (k2[0] + 2.0 * k2[1] + 2.0 * k2[2] + k2[3]) / 6.0;

	// Calulate total energy
	E = (mass / 2.0) * pow(W * L, 2) + mass * g * L * (1.0 - cos(theta));

}

/*
Runge-Kutta Formula with adaptive step size
Need the support of RKF() function
 */
void RKF_adaptive(){
	// Settings
	// ori_W        -> Original W
	// ori_theta    -> Original theta
	// yh_W         -> Record W with step h
	// yh_theta     -> Record theta with step h
	// delta        -> Indicator of truncation error
	// flag         -> 1 : delta < delta0;   0: delta > delta0 -> try again 
	GLdouble ori_W, ori_theta;
	GLdouble yh_W, yh_theta;
	GLdouble delta = 0.0;
	int flag;

	// Store original W and theta, and set flag = 0
	ori_W = W;
	ori_theta = theta;
	flag = 0;

	while(flag == 0){
		// Calculate y_h, step size = h
		RKF(dt0);
		yh_W = W;
		yh_theta = theta;

		// Calculate y_h/2, step size = h/2 -> walk 2 steps
		W = ori_W;
		theta = ori_theta;
		RKF(dt0 / 2);
		RKF(dt0 / 2);

		// calculate delta = | y2 - y1 |
		delta = fabs(theta - yh_theta) + fabs(W - yh_W);

		if(delta > delta0){
			dt0 = 0.9 * dt0 * pow(fabs(delta0 / delta), 0.25);
			flag = 0;
			
			// Recalculate y_h and y_h/2
			// Reset W and theta to origin
			W = ori_W;
			theta = ori_theta;

			continue;
		}
		else{
			flag = 1;
			dt0 = 0.9 * dt0 * pow(fabs(delta0 / delta), 0.2);
		}
	}

	// Use the Richardson Extrapolation
	W = (16 * W - yh_W) / 15.0;
	theta = (16 * theta - yh_theta) / 15.0;
	E = (mass / 2.0) * pow(W * L, 2) + mass * g * L * (1.0 - cos(theta));

}


/* Handler for window-repaint event. Call back when the window first appears and
whenever the window needs to be re-painted. */
void display(){
	
	glClearColor(1.0, 1.0, 1.0, 1.0);                       // white background
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);     // clear the color buffer
		
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, 10.0, 0, 0, 0, 0, 1, 0);                // 3d view direction
		
	glTranslatef(0.0, 0.0, 0.0);                            // translate the central point of the system
	glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex2f(-10, 0);
	glVertex2f(10, 0);
	glEnd();
		
	double x1 = L * cos(theta), y1 = L * sin(theta);	    // bob's position
	glRotatef(-90.0, 0.0, 0.0, 1.0);	                    // rotate the system clockwise by 90 degree
	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0); 
	glVertex2f(0, 0);
	glColor3f(1.0, 0.0, 0.0); 
	glVertex2f(x1, y1);
	glEnd();
	
	double r = 0.5;		                                    // bob's radius
	glBegin(GL_TRIANGLE_FAN);
	glColor3f(0.0, 0.0, 1.0);
	glVertex2f(x1, y1);
	for(double angle = 0.0; angle <= 2 * M_PI; angle = angle+0.01) {
		glColor3f(0.0, 0.0, 1.0);
		glVertex2f(x1 + sin(angle) * r, y1 + cos(angle) * r);
	}
	glEnd();
		
	glutSwapBuffers();			                  			// Double buffered - swap the front and back buffers
	
	RKF_adaptive();
	// RKF(0.05);
}

/* Handler for window re-size event. Called back when the window first appears and
	 whenever the window is re-sized with its new width and height */
void WindowSize(int w, int h){
	glViewport(0, 0, w, h);           // initial draw position and the size
	glMatrixMode(GL_PROJECTION); 	    // project 3d into 2d
	glLoadIdentity();		    // clean up the matrix data and load it to 1
	glOrtho(-10,10,-10,10,-10,10);    // matrix box of x,y,z (-10,10)
	glMatrixMode(GL_MODELVIEW); 	    // load the matrix model of rotation, scaling and translation
	glLoadIdentity();
}

/* Called back when there is no other event to be handled */
void idle(){
	glutPostRedisplay();              // post a re-paint request to activate display()
}


int main(int argc, char** argv){

	glutInit(&argc, argv);                            // initialize GLUT
	printf("Please input the initial angle (in unit of degree): ");
	scanf("%lf", &theta);
	theta = theta/180.0 * M_PI;
	printf("Please input the initial angular velocity (in unit of degree/sec): ");
	scanf("%lf", &W);
	W = W/180.0 * M_PI;

	// Calculate the initial total energy
	E0 = (mass / 2.0) * pow(W * L, 2) + mass * g * L * (1.0 - cos(theta)); // Energy of the physical pendulum
	// W = sqrt(2.0/mass/L/L*(E-mass*g*L*(1.0-cos(theta))));
	
	// GLUT
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);      // enable double buffered mode
	glutInitWindowSize(500,500);         	          // set the window's initial width & height
	glutInitWindowPosition(300,100);         	      // position the window's initial top-left corner
	glutCreateWindow("Pendulum Animation");           // create window with the given title
	glutReshapeFunc(WindowSize);                      // register callback handler for window re-size event
	glutDisplayFunc(display);                         // register callback handler for window re-paint event
	glutIdleFunc(idle);                               // register callback handler if no other event
	glutMainLoop();                                   // enter the infinite event-processing loop
 
	return 0;
} 
