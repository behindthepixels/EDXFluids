#pragma once

#include "../Base/Fluid.h"

namespace EDX
{
	namespace FluidSim
	{
		template<uint Dimension>
		class SmokeSolver : public FluidSolver<Dimension>
		{
		private:
			// Density
			DimensionalArray<Dimension, float> mDensity;
			DimensionalArray<Dimension, float> mTempDensity;
			// Temperature
			DimensionalArray<Dimension, float> mTemperature;
			DimensionalArray<Dimension, float> mTempTemperature;

			DimensionalArray<Dimension, float> mBuoyancy;

			// Smoke parameters
			float mfDiffuseRate;
			float mfHeatDiffuseRate;
			float mfTargetHeat;

			float mfBuoyDen;
			float mfBuoyHeat;


		public:
			void Initialize(const Vec<Dimension, uint>& vDim, const Vec<Dimension, float> vFluidExtent);
			void Advance(const float fDt);

			void AddSources(const float fDt);

			void WriteToFile() const;

			const float* GetDensity() { return mDensity.Data(); }
		};
	}
}