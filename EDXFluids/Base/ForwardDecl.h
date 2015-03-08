#pragma once

namespace EDX
{
	namespace FluidSim
	{
		enum Component { X = 0, Y, Z, N };
		enum CellType { Fluid = 1, Air = 1 << 2, Solid = 1 << 3, NearBoundary = 1 << 4 };
	}
}