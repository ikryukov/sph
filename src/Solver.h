#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
//#include <omp.h>
#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include "Particle.h"
#include "Grid.h"

inline float maxf(const float a, const float b)
{
	return (a>b)?a:b;
}

inline float minf(const float a, const float b)
{
	return (a<b)?a:b;
}

struct SimulationConstants
{
	float g_fSmoothlen;
	float g_fDensityCoef;

	float g_fGradPressureCoef;
	float g_fLapViscosityCoef;

	float g_fPressureStiffness;
	float g_fRestDensity;
	float g_fInvRestDensity;

	int g_iNumParticles;
	float g_fTimeStep;
	float g_fWallStiffness;
	float g_fInitialParticleSpacing;
	glm::vec4 g_vGravity;
	glm::vec4 g_vGridDim;
	glm::vec3 g_vPlanes[4];
};

class Solver
{
private:
	//std::vector<Particle> particles;
	Grid grid;
	SimulationConstants simConsts;

	//--------------------------------------------------------------------------------------
	// Density Calculation
	//--------------------------------------------------------------------------------------
	float CalculateDensity(float r_sq)
	{
		const float h_sq = simConsts.g_fSmoothlen * simConsts.g_fSmoothlen;
		// Implements this equation:
		// W_poly6(r, h) = 315 / (64 * pi * h^9) * (h^2 - r^2)^3
		// g_fDensityCoef = fParticleMass * 315.0f / (64.0f * PI * fSmoothlen^9)
		return simConsts.g_fDensityCoef * (h_sq - r_sq) * (h_sq - r_sq) * (h_sq - r_sq);
	}

	//--------------------------------------------------------------------------------------
	// Simple N^2 Algorithm
	//--------------------------------------------------------------------------------------
	void Density()
	{
		const float h_sq = simConsts.g_fSmoothlen * simConsts.g_fSmoothlen;
		#pragma omp parallel for
		for (int y = 0; y < grid.resolution.y; ++y)
		{
			for (int x = 0; x < grid.resolution.x; ++x)
			{
				std::vector<Particle>& particles = grid.getParticles(x, y);
				for (size_t k = 0; k < particles.size(); ++k) 
				{		
					glm::vec2 P_position = particles[k].position;
					float density = 0;
					const int dx[9] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
					const int dy[9] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
					// Calculate the density based on all neighbors
					for (int N_ID = 0 ; N_ID < 9 ; ++N_ID)
					{
						std::vector<Particle>& cellNeighbor = grid.getParticles(x + dx[N_ID], y + dy[N_ID]);
						for (size_t PN_ID = 0; PN_ID < cellNeighbor.size(); ++PN_ID)
						{
							//float2 N_position = cellNeighbor[PN_ID].position;
							//float2 diff = N_position - P_position;
							glm::vec2 diff = cellNeighbor[PN_ID].position - P_position;
                            float r_sq = glm::dot(diff, diff);
							if (r_sq < h_sq)
							{
								density += CalculateDensity(r_sq);
							}
						}
					}
					particles[k].density = density;
				}
			}
		}
	}

	//--------------------------------------------------------------------------------------
	// Force Calculation
	//--------------------------------------------------------------------------------------

	float CalculatePressure(const float density) const
	{
		// Implements this equation:
		// Pressure = B * ((rho / rho_0)^y  - 1)
		return simConsts.g_fPressureStiffness * maxf(pow(density * simConsts.g_fInvRestDensity, 3) - 1.0f, 0.0f);
	}

	__inline glm::vec2 CalculateGradPressure(const float r, const float P_pressure, const float N_pressure, const float N_density, const glm::vec2& diff) const
	{
		//const float h = simConsts.g_fSmoothlen;
		//float avg_pressure = 0.5f * (N_pressure + P_pressure);
		// Implements this equation:
		// W_spkiey(r, h) = 15 / (pi * h^6) * (h - r)^3
		// GRAD( W_spikey(r, h) ) = -45 / (pi * h^6) * (h - r)^2
		// g_fGradPressureCoef = fParticleMass * -45.0f / (PI * fSmoothlen^6)
		//return simConsts.g_fGradPressureCoef * avg_pressure * (h - r) * (h - r) / (N_density * r) * (diff);
		return 0.5f * simConsts.g_fGradPressureCoef * (N_pressure + P_pressure) * (simConsts.g_fSmoothlen - r) * (simConsts.g_fSmoothlen - r) / (N_density * r) * (diff);
	}

	glm::vec2 CalculateLapVelocity(const float r, const glm::vec2& P_velocity, const glm::vec2& N_velocity, const float N_density)
	{
		const float h = simConsts.g_fSmoothlen;
		//const float2 vel_diff = N_velocity - P_velocity;
		// Implements this equation:
		// W_viscosity(r, h) = 15 / (2 * pi * h^3) * (-r^3 / (2 * h^3) + r^2 / h^2 + h / (2 * r) - 1)
		// LAPLACIAN( W_viscosity(r, h) ) = 45 / (pi * h^6) * (h - r)
		// g_fLapViscosityCoef = fParticleMass * fViscosity * 45.0f / (PI * fSmoothlen^6)
		return simConsts.g_fLapViscosityCoef * (h - r) / N_density * (N_velocity - P_velocity);
	}


	//--------------------------------------------------------------------------------------
	// Simple N^2 Algorithm
	//--------------------------------------------------------------------------------------
	void Force()
	{
		const float h_sq = simConsts.g_fSmoothlen * simConsts.g_fSmoothlen;
		#pragma omp parallel for
		for (int y = 0; y < grid.resolution.y; ++y)
		{
			for (int x = 0; x < grid.resolution.x; ++x)
			{
				std::vector<Particle>& particles = grid.getParticles(x, y);
				size_t particlesCount = particles.size();
				for (size_t k = 0; k < particlesCount; ++k) 
				{	
					const glm::vec2& P_position = particles[k].position;
					const glm::vec2& P_velocity = particles[k].velocity;
					const float& P_density = particles[k].density;
					float P_pressure = CalculatePressure(P_density);

					glm::vec2 acceleration(0.0f);

					// Calculate the acceleration based on all neighbors
					const int dx[9] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
					const int dy[9] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
					// Calculate the density based on all neighbors
					for (int N_ID = 0 ; N_ID < 9 ; ++N_ID)
					{
						std::vector<Particle>& cellNeighbor = grid.getParticles(x + dx[N_ID], y + dy[N_ID]);
						size_t cellSize = cellNeighbor.size();
						for (size_t PN_ID = 0; PN_ID < cellSize; ++PN_ID)
						{
							const glm::vec2& N_position = cellNeighbor[PN_ID].position;

							const glm::vec2& diff = N_position - P_position;
                            const float r_sq = glm::dot(diff, diff);
							if (r_sq < h_sq && r_sq > 1e-5)
							{
								const glm::vec2& N_velocity = cellNeighbor[PN_ID].velocity;
								const float N_density = cellNeighbor[PN_ID].density;
								const float N_pressure = CalculatePressure(N_density);
								const float r = sqrt(r_sq);

								// Pressure Term
								acceleration += CalculateGradPressure(r, P_pressure, N_pressure, N_density, diff);

								// Viscosity Term
								acceleration += CalculateLapVelocity(r, P_velocity, N_velocity, N_density);
							}
						}

						particles[k].acceleration = acceleration / P_density;
					}
				}
			}
		}
	}

	//--------------------------------------------------------------------------------------
	// Integration
	//--------------------------------------------------------------------------------------
	void Integrate()
	{
		//#pragma omp parallel for
		for (int y = 0; y < grid.resolution.y; ++y)
		{
			for (int x = 0; x < grid.resolution.x; ++x)
			{
				std::vector<Particle>& particles = grid.getParticles(x, y);
				size_t particlesCount = particles.size();
				for (size_t k = 0; k < particlesCount; ++k) 
				{	

					glm::vec2& position = particles[k].position;
					glm::vec2& velocity = particles[k].velocity;
					glm::vec2& acceleration = particles[k].acceleration;
					
					// Apply the forces from the map walls
					for (size_t i = 0 ; i < 4 ; ++i)
					{
						glm::vec3 tmp = glm::vec3(position, 1);
                        const float dist = glm::dot(tmp, simConsts.g_vPlanes[i]);
						acceleration += minf(dist, 0.0f) * -simConsts.g_fWallStiffness * simConsts.g_vPlanes[i].xy();
					}
					//float2 diff = float2(0.5f, 0.05f) - position;
					//float dist = dot(diff, diff) - 0.025f;
					//acceleration += std::min(dist, 0.0f) * simConsts.g_fWallStiffness * 100.0f * (diff);

					// Apply gravity
					acceleration += simConsts.g_vGravity.xy();

					// Integrate
					velocity += simConsts.g_fTimeStep * acceleration;
					position += simConsts.g_fTimeStep * velocity;

					//float d = dot(float3(position, 1), simConsts.g_vPlanes[2]);
					//if (d < 0) 
					//{
					//	position.x = 0.0f;
					//	//acceleration.y = 0.0f;
					//	//velocity = float2(0.0f, 0.0f);
					//}

					// Update
					//particles[k].position = position;
					//particles[k].velocity = velocity;
				}
			}
		}
		grid.rebuildGrid();
	}

public:
	void SimulationStep()
	{
		Density();
		Force();
		Integrate();
	}
	void setConstants(SimulationConstants cbuff)
	{
		simConsts = cbuff;
	}
	void Init()
	{
		std::vector<Particle> particles;
		particles.reserve(simConsts.g_iNumParticles);
		particles.resize(simConsts.g_iNumParticles);

		const int iStartingWidth = (int)sqrt( (float)simConsts.g_iNumParticles );
		for ( int i = 0 ; i < simConsts.g_iNumParticles ; i++ )
		{
			// Arrange the particles in a nice square
			int x = i % iStartingWidth;
			int y = i / iStartingWidth;
			particles[ i ].position = glm::vec2( simConsts.g_fInitialParticleSpacing * (float)x + 0.0045f, simConsts.g_fInitialParticleSpacing * (float)y );
		}
		grid.initByResolution(glm::ivec2(128, 128), glm::vec2(1.0, 1.0));
		grid.setParticles(&particles[0], particles.size());
	}
	SimulationConstants& getSimulationConstants()
	{
		return simConsts;
	}
	Grid& getGrid()
	{
		return grid;
	}
	Solver(void);
	~Solver(void);
};

