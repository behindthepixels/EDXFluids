#include "Smoke.h"
#include "Math/Vector.h"

#include <fstream>

namespace EDX
{
	namespace FluidSim
	{
		template<uint Dimension>
		void SmokeSolver<Dimension>::Initialize(const Vec<Dimension, uint>& fluidDim, const Vec<Dimension, float> vFluidExtent)
		{
			mbUseFLIP = false;
			mbReseedParticles = false;
			mbUseVorticityConfinement = false;
			mbRebuildPossionPerFrame = false;

			FluidSolver<Dimension>::Initialize(fluidDim, vFluidExtent);

			mDensity.Init(fluidDim);
			mTempDensity.Init(fluidDim);
			mTemperature.Init(fluidDim);
			mTempTemperature.Init(fluidDim);
			mBuoyancy.Init(fluidDim);

			mfDiffuseRate = 0.0f;
			mfTargetHeat = 8.4f;
			mfHeatDiffuseRate = 0.00005f;

			mfBuoyDen = 0.0625f * 0.5f;
			mfBuoyHeat = 0.025f;
		}

		template<uint Dimension>
		void SmokeSolver<Dimension>::Advance(float fDt)
		{
			FluidSolver<Dimension>::Advance(fDt);

			AdvectValueMacCormack(fDt, mDensity, mTempDensity);
			DiffuseValueGaussSedel(fDt, mfDiffuseRate, mDensity, mTempDensity);

			AdvectValueMacCormack(fDt, mTemperature, mTempTemperature);
			DiffuseValueGaussSedel(fDt, mfHeatDiffuseRate, mTemperature, mTempTemperature);

			if(Dimension == 3)
				WriteToFile();
		}

		template<uint Dimension>
		void SmokeSolver<Dimension>::AddSources(const float fDt)
		{
			// Add sources
			parallel_for(0, (int)mDensity.LinearSize(), [&](int i)
			{
				Vec<Dimension, float> vPos = mDensity.Index(i);
				if(vPos[1] >= 0 && vPos[1] <= mvDim[1] / 20.0f)
				{
					Vec<Dimension, float> vCen = 0.5f * Vec<Dimension, float>(mvDim);
					vCen[1] = vPos[1];

					float fDist = Math::Length(vPos - vCen);
					if(fDist <= mvDim[1] / 15.0f)
					{
						mTemperature[i] += (1.0f - Math::Exp(-mfTargetHeat * fDt)) * (mfTargetHeat - mTemperature[i]);
						mDensity[i] = 1.0f;
					}
					else
					{
						mTemperature[i] += (1.0f - Math::Exp(0.0f * fDt)) * (mfTargetHeat - mTemperature[i]);
					}
				}
				else
				{
					mTemperature[i] += (1.0f - Math::Exp(0.0f * fDt)) * (mfTargetHeat - mTemperature[i]);
				}
			});
// 			for(auto i = mvDim[0]/2 - mvDim[0]/15; i < mvDim[0]/2 + mvDim[0]/15; i++)
// 			{
// 				for(auto j = mvDim[1]/20; j < mvDim[1]/10; j++)
// 				{
// 					mDensity[Vec<Dimension, uint>(i, j)] = 1.0f;
// 				}
// 			}

// 			parallel_for(0, (int)mTemperature.LinearSize(), [&](int i)
// 			{
// 				Vec<Dimension, uint> vIdx = mTemperature.Index(i);
// 				if(vIdx.x >= mvDim[0]/2 - mvDim[0]/15 && vIdx.x < mvDim[0]/2 + mvDim[0]/15 &&
// 					vIdx.y >= mvDim[1]/20 && vIdx.y < mvDim[1]/10)
// 					mTemperature[i] += (1.0f - Math::Exp(-mfTargetHeat * fDt)) * (mfTargetHeat - mTemperature[i]);
// 				else
// 					mTemperature[i] += (1.0f - Math::Exp(0.0f * fDt)) * (mfTargetHeat - mTemperature[i]);
// 			});

			parallel_for(0, (int)mBuoyancy.LinearSize(), [&](int i)
			{
				mBuoyancy[i] = mfBuoyDen * mDensity[i] - mfBuoyHeat * (mTemperature[i] - 0.0f);
			});

			parallel_for(0, (int)mVelocity[1].LinearSize(), [&](int i)
			{
				Vec<Dimension, float> vOffset;
				vOffset[1] = -0.5f;
				Vec<Dimension, uint> vIdx = mVelocity[1].Index(i);
				if(mFluidVolume[1][i] == 0.0f)
					return;

				Vec<Dimension, uint> vUpIdx = vIdx;
				vUpIdx[1]--;
				mVelocity[1][i] -= 0.5f * (mBuoyancy[vUpIdx] + mBuoyancy[vIdx]) * fDt * mfInvDx;
			});

			//mVelocity[1][Vec<Dimension, uint>(mvDim[0]/2, mvDim[1]*0.07f)] = 2.0f * mfInvDx;
		}

		template<uint Dimension>
		void SmokeSolver<Dimension>::WriteToFile() const
		{
			char strFileName[260];
			sprintf_s(strFileName, 260, "../Output/FluidData%i.txt", miFrameIdx);

			std::ofstream outFile;
			outFile.open(strFileName);
			assert(outFile.is_open());

			for(int i = 0; i < mDensity.LinearSize(); i++)
				outFile << mDensity[i] << " ";

			outFile.close();
		}


		template class SmokeSolver<2>;
		template class SmokeSolver<3>;
	}
}