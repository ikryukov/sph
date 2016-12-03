#include <GL/glew.h>
#include <GLUT/GLUT.h>
#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <math.h>
#include "Solver.h"
#include "Grid.h"
#ifdef _WINDOWS
#include <Windows.h>
#endif
#define MAX_FRAMES 1 << 10

int frameNumber = 0;
Solver solver;

static glm::vec3 Rainbow[5] = {
    glm::vec3(1, 0, 0), // red
    glm::vec3(1, 1, 0), // orange
    glm::vec3(0, 1, 0), // green
    glm::vec3(0, 1, 1), // teal
    glm::vec3(0, 0, 1), // blue
};

float saturate(float n)
{
	float res = n;
	if (res > 1) res = 1.0f;
	if (res < 0) res = 0.0f;
	return res;
}

glm::vec3 VisualizeNumber(float n)
{
	float a = floorf(n * 4.0f);
	float b = ceilf(n * 4.0f);
	float c = saturate(b - a);
    glm::vec3 t = (Rainbow[ (int) b ] - Rainbow[ (int) a]);
    t *= c;
	return Rainbow[ (int) a] + t;
}

glm::vec3 VisualizeNumber(float n, float lower, float upper)
{
    return VisualizeNumber( saturate( (n - lower) / (upper - lower) ) );
}


void renderScene(void) 
{		
	glLoadIdentity();
	glTranslatef(-0.7, -0.7, 0.0);
	glScalef(2.0, 2.0, 1.0);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glColor3f( 0.552f, 0.713f, 0.803f);
	glBegin(GL_QUADS);
	glVertex3f(-1.0, -1.0, 0.0);
	glVertex3f(-1.0, 1.0, 0.0);
	glVertex3f(1.0, 1.0, 0.0);
	glVertex3f(1.0, -1.0, 0.0);
	glEnd();

	glColor3f( 1.0f, 1.0f, 1.0f);
	glBegin(GL_LINES);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.7, 0.0, 0.0);
	glVertex3f(0.7, 0.0, 0.0);
	glVertex3f(0.7, 0.7, 0.0);
	glVertex3f(0.7, 0.7, 0.0);
	glVertex3f(0.0, 0.7, 0.0);
	glVertex3f(0.0, 0.7, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glEnd();

	//glColor3f( 0.852f, 0.13f, 0.03f);
	glBegin(GL_POINTS);
	glPointSize(2.0);
	Grid grid = solver.getGrid();
	for (int y = 0; y < grid.resolution.y; ++y)
	{
		for (int x = 0; x < grid.resolution.x; ++x)
		{
			std::vector<Particle>& particles = grid.getParticles(x, y);
			for (size_t k = 0; k < particles.size(); ++k) 
			{
				glm::vec3 color = VisualizeNumber(particles[k].density, 1000.0f, 2000.0f);
				glColor3f(color.x, color.y, color.z);
				glVertex3f(particles[k].position.x, particles[k].position.y, 0);
			}
		}
	}
	glEnd();
	glutSwapBuffers();
	float simulation_time = solver.getSimulationConstants().g_fTimeStep * frameNumber;
	printf("Frame: %d \t Simulation time: %.4f seconds \n", frameNumber, simulation_time);
#ifdef _WINDOWS
	long long _start_tick = GetTickCount();	
#endif
	solver.SimulationStep();
#ifdef _WINDOWS
	long long _finish_tick = GetTickCount();
	long long _elapsed_tick = _finish_tick - _start_tick;
	printf("Done in %d milliseconds\n", _elapsed_tick);
#endif
	if (frameNumber > MAX_FRAMES) 
	{
		exit(0);
	}
	++frameNumber;
}

void processNormalKeys(unsigned char key, int x, int y) 
{
	if (key == 27)
		exit(0);
}

int main(int argc, char** argv)
{
	float g_fInitialParticleSpacing = 0.0045f;
	float g_fSmoothlen = 0.012f;
	float g_fPressureStiffness = 200.0f;
	float g_fRestDensity = 1000.0f;
	float g_fParticleMass = 0.0002f;
	float g_fViscosity = 0.1f;
	float g_fMaxAllowableTimeStep = 0.005f;
	float g_fParticleRenderSize = 0.003f;
	float PI = 3.1415926;
	float g_fWallStiffness = 3000.0f;

	float g_fMapHeight = 0.7f;
	//float g_fMapWidth = (4.0f / 3.0f) * g_fMapHeight;
	float g_fMapWidth = 0.7f;
	glm::vec3 g_vPlanes[4] = {
		glm::vec3(1, 0, 0),
		glm::vec3(0, 1, 0),
		glm::vec3(-1, 0, g_fMapWidth),
		glm::vec3(0, -1, g_fMapHeight)
	};

	const int NUM_PARTICLES_1K = 1 * 1024;
	const int NUM_PARTICLES_2K = 2 * 1024;
	const int NUM_PARTICLES_4K = 4 * 1024;
	const int NUM_PARTICLES_8K = 8 * 1024;
	const int NUM_PARTICLES_16K = 16 * 1024;
	const int NUM_PARTICLES_32K = 32 * 1024;
	const int NUM_PARTICLES_64K = 64 * 1024;
	int g_iNumParticles = NUM_PARTICLES_4K;

	// Gravity Directions
	const glm::vec2 GRAVITY_DOWN(0, -0.5f);
	const glm::vec2 GRAVITY_UP(0, 0.5f);
	const glm::vec2 GRAVITY_LEFT(-0.5f, 0);
	const glm::vec2 GRAVITY_RIGHT(0.5f, 0);
	glm::vec2 g_vGravity = GRAVITY_DOWN;

	SimulationConstants cbuff;
	cbuff.g_fDensityCoef = g_fParticleMass * 315.0f / (64.0f * PI * pow(g_fSmoothlen, 9));
	cbuff.g_fGradPressureCoef = g_fParticleMass * -45.0f / (PI * pow(g_fSmoothlen, 6));
	cbuff.g_fInitialParticleSpacing = g_fInitialParticleSpacing;
	cbuff.g_fLapViscosityCoef = g_fParticleMass * g_fViscosity * 45.0f / (PI * pow(g_fSmoothlen, 6));
	cbuff.g_fPressureStiffness = g_fPressureStiffness;
	cbuff.g_fRestDensity = g_fRestDensity;
	cbuff.g_fInvRestDensity = 1.0f / g_fRestDensity;
	cbuff.g_fSmoothlen = g_fSmoothlen;
	cbuff.g_fTimeStep = g_fMaxAllowableTimeStep;
	cbuff.g_fWallStiffness = g_fWallStiffness;
	cbuff.g_iNumParticles = g_iNumParticles;
	cbuff.g_vGravity = glm::vec4(g_vGravity.x, g_vGravity.y, 0, 0);
	cbuff.g_vGridDim.x = 1.0f / g_fSmoothlen;
	cbuff.g_vGridDim.y = 1.0f / g_fSmoothlen;
	cbuff.g_vGridDim.z = 0;
	cbuff.g_vGridDim.w = 0;
	cbuff.g_vPlanes[0] = g_vPlanes[0];
	cbuff.g_vPlanes[1] = g_vPlanes[1];
	cbuff.g_vPlanes[2] = g_vPlanes[2];
	cbuff.g_vPlanes[3] = g_vPlanes[3];
	
	solver.setConstants(cbuff);
	solver.Init();

	glutInit(&argc, argv);
	glutInitWindowPosition(50, 50);
	glutInitWindowSize(800, 800);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("SPH");
	glutDisplayFunc(renderScene);
	glutIdleFunc(renderScene);
	glutKeyboardFunc(processNormalKeys);	
	glutMainLoop();
	return 0;
}
