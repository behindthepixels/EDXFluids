#pragma once

#include "../Base/Fluid.h"
#include "LevelSet.h"

namespace EDX
{
	namespace FluidSim
	{
		template<uint Dimension>
		class LiquidSolver : public FluidSolver<Dimension>
		{
		public:
			void Initialize(const Vec<Dimension, uint>& vDim, const Vec<Dimension, float> vFluidExtent);
			void Advance(const float fDt);
			void Step(const float fDt);
			
		private:
			void Extrapolate();
			void AddSources(const float fDt);

			//void InitMarkers();
			void UpdateMarkers();
			
			void WriteToFile() const;
		};
	}
}