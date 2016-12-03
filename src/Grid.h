#pragma once
#include <cstdlib>
#include <cmath>
#include <vector>
#include <glm/glm.hpp>
#include "Particle.h"

struct Bucket {
	std::vector<Particle> p;
	void addParticle(Particle& particle)
	{
		p.push_back(particle);
	}
	Bucket()
	{
//		p.reserve(16);
//		p.resize(16);
	}
};

class Grid
{
public:
    glm::vec2 min; //min point in world space
	glm::vec2 max; //max point in world space
	glm::vec2 step; // step by x, y in world space
	glm::ivec2 resolution;
	std::vector<Particle> empty;
	Bucket* cells;
	// ^y
	// |  local coordinates, row-major in memory
	// +----> x
	// 0
public:
	Grid(void);
	void initByStep(glm::vec2 steps, glm::vec2 area)
	{
		step = steps;
		resolution.x = (int) ceilf(area.x / steps.x);
		resolution.y = (int) ceilf(area.y / steps.y);
		cells = new Bucket[resolution.y * resolution.x];
	}
	void initByResolution(glm::ivec2 resolutions, glm::vec2 area)
	{
		resolution = resolutions;
		step.x = area.x / resolutions.x;
		step.y = area.y / resolutions.y;
		cells = new Bucket[resolution.y * resolution.x];
	}
	int getCellIndex(glm::vec2 pos)
	{
		return abs((int) (pos.y / step.y)) * resolution.x + abs((int) (pos.x / step.x));
	}
	void setParticles(Particle* p, int size)
	{
		for (int i = 0; i < size; ++i)
		{
			int cellIdx = getCellIndex(p[i].position);
			cells[cellIdx].addParticle(p[i]);
		}
	}
	void rebuildGrid()
	{
		for (int i = 0; i < resolution.y; ++i)
		{
			for (int j = 0; j < resolution.x; ++j)
			{
				int refCellIdx = i * resolution.x + j;
				for (std::vector<Particle>::iterator k = cells[refCellIdx].p.begin() ; k != cells[refCellIdx].p.end(); )
				{
					int cellIdx = getCellIndex(k->position);
					if (cellIdx != refCellIdx)
					{
						cells[cellIdx].addParticle(*k);
						k = cells[refCellIdx].p.erase(k);
					} 
					else
					{
						++k;
					}
				}
			}
		}
	}
	std::vector<Particle>& getParticles(int cellX, int cellY)
	{
		if (cellX < 0 || cellY < 0 || cellX >= resolution.x || cellY >= resolution.y)
		{
			return empty;
		}
		return cells[cellY * resolution.x + cellX].p;
	}
	~Grid(void);
};

