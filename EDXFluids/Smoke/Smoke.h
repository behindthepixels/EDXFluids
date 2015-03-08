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
			Array<Dimension, float> mDensity;
			Array<Dimension, float> mTempDensity;
			// Temperature
			Array<Dimension, float> mTemperature;
			Array<Dimension, float> mTempTemperature;

			Array<Dimension, float> mBuoyancy;

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