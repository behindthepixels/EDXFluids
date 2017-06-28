#pragma once

#include "EDXPrerequisites.h"
#include "../Base/ForwardDecl.h"
#include "Containers/DimensionalArray.h"
#include "Math/Vector.h"

namespace EDX
{
	namespace FluidSim
	{
		template<uint Dimension>
		class LevelSet
		{
		protected:
			DimensionalArray<Dimension, float> mPhi;
			DimensionalArray<Dimension, float> mPhiPrev;
			DimensionalArray<Dimension, Vec<Dimension, int>> mClosestSurfIdx;

		public:
			void Initialize(const Vec<Dimension, int>& vDim);
			template<typename PhiFunc>
			void InitSDF(PhiFunc func)
			{
				for(auto i = 0; i < mPhi.LinearSize(); i++)
				{
					Vec<Dimension, float> vPos = mPhi.Index(i);
					mPhi[i] = func(vPos);
				}
			}
			void ReInitSDF(const DimensionalArray<Dimension, CellType>& markers, const CellType type, const float fNarrowBand = 1e10f);
			const DimensionalArray<Dimension, float>& GetPhi() const { return mPhi; }
			Vec<Dimension, int> GetClosestSurfIdx(const Vec<Dimension, int>& vIdx) const { return mClosestSurfIdx[vIdx]; }
		};
	}
}