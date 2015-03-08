#include "Liquid.h"
#include "Math/Vector.h"

namespace EDX
{
	namespace FluidSim
	{
		template<uint Dimension>
		void LiquidSolver<Dimension>::Initialize(const Vec<Dimension, uint>& vDim, const Vec<Dimension, float> vFluidExtent)
		{
			mbUseFLIP = true;
			mbReseedParticles = false;
			mbUseVorticityConfinement = false;
			mbRebuildPossionPerFrame = true;

			FluidSolver<Dimension>::Initialize(vDim, vFluidExtent);
		}
		
		template<uint Dimension>
		void LiquidSolver<Dimension>::Advance(const float fDt)
		{
			FluidSolver<Dimension>::Advance(fDt);
			
			if(Dimension == 3)
				WriteToFile();
		}
		
		template<uint Dimension>
		void LiquidSolver<Dimension>::Step(const float fDt)
		{
			FluidSolver<Dimension>::Step(fDt);
			Extrapolate();
		}

		template<uint Dimension>
		void LiquidSolver<Dimension>::Extrapolate()
		{
			for(auto d = 0; d < Dimension; d++)
			{
				parallel_for(0, (int)mVelocity[d].LinearSize(), [&](int i)
				{
					Vec<Dimension, int> vIdx1 = mVelocity[d].Index(i);
					Vec<Dimension, int> vIdx2 = vIdx1 - Vec<Dimension, int>::UNIT[d];
					vIdx1[d] = Math::Clamp(vIdx1[d], 0, mvDim[d] - 1);
					vIdx2[d] = Math::Clamp(vIdx2[d], 0, mvDim[d] - 1);
					if((mMarkers[vIdx1] & Fluid) && (mMarkers[vIdx2] & Fluid))
						return;

					const Vec<Dimension, int> vSurfIdx1 = mLevelSet.GetClosestSurfIdx(vIdx1);
					const Vec<Dimension, int> vSurfIdx2 = mLevelSet.GetClosestSurfIdx(vIdx2);
					
					if(vSurfIdx1 == 0)
						return;

					const float fVel1 = GetValueFace(vSurfIdx1, mVelocity, Component(d));
					const float fVel2 = GetValueFace(vSurfIdx2, mVelocity, Component(d));

					mVelocity[d][i] = 0.5f * (fVel1 + fVel2);
				});
			}
		}

		template<uint Dimension>
		void LiquidSolver<Dimension>::AddSources(const float fDt)
		{
			parallel_for(0, (int)mVelocity[1].LinearSize(), [&](int i)
			{
				Vec<Dimension, float> vOffset;
				vOffset[1] = -0.5f;
				Vec<Dimension, uint> vIdx = mVelocity[1].Index(i);
				if(GetValueFace(Vec<Dimension, float>(vIdx) + vOffset, mVolume) == Vec<Dimension, float>::ZERO)
					return;

				mVelocity[1][i] -= 0.3f * fDt * mfInvDx;;
			});

		}

		template<uint Dimension>
		void LiquidSolver<Dimension>::UpdateMarkers()
		{
			parallel_for(0, (int)mMarkers.LinearSize(), [&](int i)
			{
				if(!(mMarkers[i] & Solid))
					mMarkers[i] = Air;
			});

			for(const auto& particle : mParticles)
			{
				const Vec<Dimension, int> vIdx = Math::RoundToInt(particle.vPos);

				if(!(mMarkers[vIdx] & Solid))
					mMarkers[vIdx] = Fluid;
			}

			mLevelSet.ReInitSDF(mMarkers, Fluid, mvDim[1] * 0.02f);

			// Calculate near boundary flags
			SetBoundaryMarkers();

			for(auto dim = 0; dim < Dimension; dim++)
			{
				parallel_for(0, (int)mFluidVolume[dim].LinearSize(), [&](int i)
				{
					Vec<Dimension, uint> vIdx = mFluidVolume[dim].Index(i);
					if(vIdx[dim] <= 1 || vIdx[dim] >= mLevelSet.GetPhi().Size(dim) - 1)
					{
						mFluidVolume[dim][vIdx] = 0.0f;
						return;
					}
					for(auto d = 0; d < Dimension; d++)
					{
						if(d == dim)
							continue;
						if(vIdx[d] <= 0 || vIdx[d] >= mLevelSet.GetPhi().Size(d) - 1)
						{
							mFluidVolume[dim][vIdx] = 0.0f;
							return;
						}
					}
					
					const Vec<Dimension, uint> vOffset = Vec<Dimension, uint>::UNIT[dim];
					mFluidVolume[dim][i] = Math::Clamp(FractionInside(mLevelSet.GetPhi()[vIdx - vOffset], mLevelSet.GetPhi()[vIdx]), 0.0f, 1.0f);
				});
			}
		}
		
		template<uint Dimension>
		void LiquidSolver<Dimension>::WriteToFile() const
		{
			char strFileName[260];
			sprintf_s(strFileName, 260, "../Output/ParticleData%i.txt", miFrameIdx);

			std::ofstream outFile;
			outFile.open(strFileName);
			assert(outFile.is_open());

			for(const auto& particle : mParticles)
				outFile << particle.vPos[0] << " " << particle.vPos[1] << " " << particle.vPos[2] << "\n";

			outFile.close();
		}

		template class LiquidSolver<2>;
		template class LiquidSolver<3>;
	}
}