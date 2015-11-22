#include "Windows/Window.h"
#include "Windows/Application.h"
#include "Graphics/EDXGui.h"
#include "Graphics/OpenGL.h"
#include "Graphics/Camera.h"

#include "Base/Fluid.h"
#include "Smoke/Smoke.h"
#include "Water/Liquid.h"
#include "Water/IsoSurface.h"
#include "Util/VectorVisualizer.h"

using namespace EDX;
using namespace EDX::FluidSim;
using namespace EDX::GUI;

LiquidSolver<3> gFluid;

int giDim = 192;
const float gfTimeStep = 1.0f / 30.0f;
bool gbSimulate = true;
int debugMode = 0;

Camera gCamera;

void OnInit(Object* pSender, EventArgs args)
{
	//IsoSurface isoSurface;
	//isoSurface.Initialize(192, 192, 192, 1);

	//std::ifstream inFile;
	//inFile.open("C:/Users/Edward Liu/Dropbox/LevelSetData66.txt");

	//const auto size = 192 * 192 * 192;
	//float* data = new float[size];

	//Array3f sphere;
	//sphere.Init(Vector3i(192, 192, 192));
	//for (auto i = 0; i < sphere.LinearSize(); i++)
	//{
	//	Vector3 pos = sphere.Index(i);
	//	float dist = Math::Distance(pos, Vector3(96, 96, 96));
	//	data[i] = dist < 48.0f ? -1e5f: 1e5f;
	//	inFile >> data[i];
	//}

	//isoSurface.UpdateSurface(data);
	//isoSurface.MarchingCube();
	//isoSurface.GenerateObjMesh();

	//SafeDeleteArray(data);

	OpenGL::InitializeOpenGLExtensions();

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	gFluid.Initialize(Vector3i(giDim, giDim, giDim), Vector3::UNIT_SCALE);
	gCamera = Camera(Vector3(96, 96, 525), Vector3(96, 96, 0), Vector3::UNIT_Y, 840, 640);

	EDXGui::Init();
}

void OnRender(Object* pSender, EventArgs args)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	if (gbSimulate)
		gFluid.Advance(gfTimeStep);
	//gVisual.AdvectRK4(gFluid.mAvgVelocity[0].Data(), gFluid.mAvgVelocity[1].Data(), 0.75f, 1.0f, gfTimeStep);
	//gVisual.Render(640.0f / float(giDim), gFluid.mVisual.Data());
	//glPointSize(1.0f);
	//ParticleVisualizer<2>::Render(gFluid.GetParticles(), 640.0f / float(giDim));

	//float fLen = 640.0f / float(giDim);

	//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glBegin(GL_QUADS);
	//for (auto i = 0; i < gFluid.mMarkers.LinearSize(); i++)
	//{
	//	auto vIdx = gFluid.mMarkers.Index(i);

	//	glColor3f(gFluid.mLevelSet.GetPhi()[vIdx], gFluid.mLevelSet.GetPhi()[vIdx], gFluid.mLevelSet.GetPhi()[vIdx]);
	//	glVertex2f(vIdx.x * fLen, vIdx.y * fLen);
	//	//glColor3f(gFluid.mLevelSet.GetPhi()[vIdx + Vector2i(1, 0)], gFluid.mLevelSet.GetPhi()[vIdx + Vector2i(1, 0)], gFluid.mLevelSet.GetPhi()[vIdx + Vector2i(1, 0)]);
	//	glVertex2f(vIdx.x * fLen + fLen, vIdx.y * fLen);
	//	//glColor3f(gFluid.mLevelSet.GetPhi()[vIdx + Vector2i(1, 1)], gFluid.mLevelSet.GetPhi()[vIdx + Vector2i(1, 1)], gFluid.mLevelSet.GetPhi()[vIdx + Vector2i(1, 1)]);
	//	glVertex2f(vIdx.x * fLen + fLen, vIdx.y * fLen + fLen);
	//	//glColor3f(gFluid.mLevelSet.GetPhi()[vIdx + Vector2i(0, 1)], gFluid.mLevelSet.GetPhi()[vIdx + Vector2i(0, 1)], gFluid.mLevelSet.GetPhi()[vIdx + Vector2i(0, 1)]);
	//	glVertex2f(vIdx.x * fLen, vIdx.y * fLen + fLen);
	//}
	//glEnd();

	glMatrixMode(GL_MODELVIEW);
	const Matrix& mView = gCamera.GetViewMatrix();
	glLoadTransposeMatrixf((float*)&mView);

	ParticleVisualizer<3>::Render(gFluid.GetParticles(), 1.0f);
}

void OnResize(Object* pSender, ResizeEventArgs args)
{
	glViewport(0, 0, args.Width, args.Height);

	// Set opengl params
	glMatrixMode(GL_PROJECTION);
	const Matrix& mProj = gCamera.GetProjMatrix();
	glLoadTransposeMatrixf((float*)&mProj);

	EDXGui::Resize(args.Width, args.Height);
}

void OnMouseEvent(Object* pSender, MouseEventArgs args)
{
	if (EDXGui::HandleMouseEvent(args))
		return;

	gCamera.HandleMouseMsg(args);
}

void OnKeyboardEvent(Object* pSender, KeyboardEventArgs args)
{
	if (EDXGui::HandleKeyboardEvent(args))
		return;

	gCamera.HandleKeyboardMsg(args);
}

void OnRelease(Object* pSender, EventArgs args)
{
	EDXGui::Release();
}


int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR cmdArgs, int cmdShow)
{
	Application::Init(hInst);

	Window* mainWindow = new GLWindow;
	mainWindow->SetMainLoop(NotifyEvent(OnRender));
	mainWindow->SetInit(NotifyEvent(OnInit));
	mainWindow->SetResize(ResizeEvent(OnResize));
	mainWindow->SetMouseHandler(MouseEvent(OnMouseEvent));
	mainWindow->SetkeyboardHandler(KeyboardEvent(OnKeyboardEvent));
	mainWindow->SetRelease(NotifyEvent(OnRelease));

	mainWindow->Create(L"EDXFluids", 840, 640);

	Application::Run(mainWindow);

	return 0;
}