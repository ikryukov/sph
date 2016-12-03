/*
 * Copyright (C) 2009 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// OpenGL ES 2.0 code
#ifdef __ANDROID__

#include <jni.h>
#include <android/log.h>

#include <GLES2/gl2.h>
#include <GLES2/gl2ext.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Solver.h"

#define  LOG_TAG    "libgl2jni"
#define  LOGI(...)  __android_log_print(ANDROID_LOG_INFO,LOG_TAG,__VA_ARGS__)
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR,LOG_TAG,__VA_ARGS__)

static void printGLString(const char *name, GLenum s) {
    const char *v = (const char *) glGetString(s);
    LOGI("GL %s = %s\n", name, v);
}

static void checkGlError(const char* op) {
    for (GLint error = glGetError(); error; error
            = glGetError()) {
        LOGI("after %s() glError (0x%x)\n", op, error);
    }
}

static const char gVertexShader[] = 
    "attribute vec4 vPosition;\n"
    "void main() {\n"
    "  gl_Position = vPosition;\n"
    "}\n";

static const char gFragmentShader[] = 
    "precision mediump float;\n"
    "void main() {\n"
    "  gl_FragColor = vec4(0.0, 1.0, 0.0, 1.0);\n"
    "}\n";

GLuint loadShader(GLenum shaderType, const char* pSource) {
    GLuint shader = glCreateShader(shaderType);
    if (shader) {
        glShaderSource(shader, 1, &pSource, NULL);
        glCompileShader(shader);
        GLint compiled = 0;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
        if (!compiled) {
            GLint infoLen = 0;
            glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLen);
            if (infoLen) {
                char* buf = (char*) malloc(infoLen);
                if (buf) {
                    glGetShaderInfoLog(shader, infoLen, NULL, buf);
                    LOGE("Could not compile shader %d:\n%s\n",
                            shaderType, buf);
                    free(buf);
                }
                glDeleteShader(shader);
                shader = 0;
            }
        }
    }
    return shader;
}

GLuint createProgram(const char* pVertexSource, const char* pFragmentSource) {
    GLuint vertexShader = loadShader(GL_VERTEX_SHADER, pVertexSource);
    if (!vertexShader) {
        return 0;
    }

    GLuint pixelShader = loadShader(GL_FRAGMENT_SHADER, pFragmentSource);
    if (!pixelShader) {
        return 0;
    }

    GLuint program = glCreateProgram();
    if (program) {
        glAttachShader(program, vertexShader);
        checkGlError("glAttachShader");
        glAttachShader(program, pixelShader);
        checkGlError("glAttachShader");
        glLinkProgram(program);
        GLint linkStatus = GL_FALSE;
        glGetProgramiv(program, GL_LINK_STATUS, &linkStatus);
        if (linkStatus != GL_TRUE) {
            GLint bufLength = 0;
            glGetProgramiv(program, GL_INFO_LOG_LENGTH, &bufLength);
            if (bufLength) {
                char* buf = (char*) malloc(bufLength);
                if (buf) {
                    glGetProgramInfoLog(program, bufLength, NULL, buf);
                    LOGE("Could not link program:\n%s\n", buf);
                    free(buf);
                }
            }
            glDeleteProgram(program);
            program = 0;
        }
    }
    return program;
}

GLuint gProgram;
GLuint gvPositionHandle;
Solver solver;

void initFluid()
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
	float3 g_vPlanes[4] = {
		float3(1, 0, 0),
		float3(0, 1, 0),
		float3(-1, 0, g_fMapWidth),
		float3(0, -1, g_fMapHeight)
	};

	const int NUM_PARTICLES_1K = 1 * 1024;
	const int NUM_PARTICLES_2K = 2 * 1024;
	const int NUM_PARTICLES_4K = 4 * 1024;
	const int NUM_PARTICLES_8K = 8 * 1024;
	const int NUM_PARTICLES_16K = 16 * 1024;
	const int NUM_PARTICLES_32K = 32 * 1024;
	const int NUM_PARTICLES_64K = 64 * 1024;
	int g_iNumParticles = NUM_PARTICLES_16K;

	// Gravity Directions
	const float2 GRAVITY_DOWN(0, -0.5f);
	const float2 GRAVITY_UP(0, 0.5f);
	const float2 GRAVITY_LEFT(-0.5f, 0);
	const float2 GRAVITY_RIGHT(0.5f, 0);
	float2 g_vGravity = GRAVITY_DOWN;

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
	cbuff.g_vGravity = float4(g_vGravity.x, g_vGravity.y, 0, 0);
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
}

bool setupGraphics(int w, int h) {
    printGLString("Version", GL_VERSION);
    printGLString("Vendor", GL_VENDOR);
    printGLString("Renderer", GL_RENDERER);
    printGLString("Extensions", GL_EXTENSIONS);

    LOGI("setupGraphics(%d, %d)", w, h);
    gProgram = createProgram(gVertexShader, gFragmentShader);
    if (!gProgram) {
        LOGE("Could not create program.");
        return false;
    }
    gvPositionHandle = glGetAttribLocation(gProgram, "vPosition");
    checkGlError("glGetAttribLocation");
    LOGI("glGetAttribLocation(\"vPosition\") = %d\n",
            gvPositionHandle);

    glViewport(0, 0, w, h);
    checkGlError("glViewport");
    initFluid();
    return true;
}

const GLfloat gTriangleVertices[] = { 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 1.0f };

void renderFrame() {
    static float grey;
    grey += 0.01f;
    if (grey > 1.0f) {
        grey = 0.0f;
    }

    glClearColor(grey, grey, grey, 1.0f);
    checkGlError("glClearColor");
    glClear( GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    checkGlError("glClear");

    glUseProgram(gProgram);
    checkGlError("glUseProgram");

    glVertexAttribPointer(gvPositionHandle, 2, GL_FLOAT, GL_FALSE, 0, gTriangleVertices);
    checkGlError("glVertexAttribPointer");
    glEnableVertexAttribArray(gvPositionHandle);
    checkGlError("glEnableVertexAttribArray");
    glDrawArrays(GL_TRIANGLES, 0, 3);
    checkGlError("glDrawArrays");
    solver.SimulationStep();
}

extern "C" {
    JNIEXPORT void JNICALL Java_ru_total_SPHLib_init(JNIEnv * env, jobject obj,  jint width, jint height);
    JNIEXPORT void JNICALL Java_ru_total_SPHLib_step(JNIEnv * env, jobject obj);
};

JNIEXPORT void JNICALL Java_ru_total_SPHLib_init(JNIEnv * env, jobject obj,  jint width, jint height)
{
    setupGraphics(width, height);
}

JNIEXPORT void JNICALL Java_ru_total_SPHLib_step(JNIEnv * env, jobject obj)
{
    renderFrame();
}
#endif // __ANDROID__
