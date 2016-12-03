#pragma once
#include <glm/glm.hpp>

class Particle
{
public:
	glm::vec2 position;
	glm::vec2 velocity;
	float density;
	glm::vec2 acceleration;
	
	//float padding;
	Particle(void);
	~Particle(void);
};

