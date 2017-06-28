#include "EDXPrerequisites.h"
#include "Math/Vector.h"
#include "Containers/DimensionalArray.h"
#include "Core/Random.h"

namespace EDX
{
	namespace Math { class FFT; }

	class Ocean
	{
	private:
		// Ocean parameters
		int miMapDimFFT;
		float mfPatchSize;
		float mfWaveAmplitude;
		float mfWindSpeed;
		float fTime;
		float mfChoppyScale;
		float mfTimeScale;
		Vector2 mvWindDir;

		float mfAmpFactor;

		// Wave frequency
		DimensionalArray<2, Vector2> mavH0;
		DimensionalArray<2, Vector2> mavHt;
		DimensionalArray<2, Vector2> mavTroppyX;
		DimensionalArray<2, Vector2> mavTroppyZ;
		Array2f mafOmega;

		Array2f mafHeightData;
		Array2f mafTroppyXData;
		Array2f mafTroppyZData;

		// RNG
		RandomGen mRNG;
		int miFrameIdx;

		static const int GRAV_CONST = 981;

	public:
		void Initialize();
		void UpdateHeightMap(const Math::FFT& fft, const float fTime);
		void GenerateObjMesh() const;

		const float* GetHeightMap() const { return mafHeightData.Data(); }
		int GetFFTDim() const { return miMapDimFFT; }

	private:
		float PhillipsSpectrum(const Vector2& vK) const;
		void InitHeightMap();
		Vector3 GetVertex(const Vector2i& vIdx) const;
	};
}