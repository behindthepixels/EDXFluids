#include "Windows/Window.h"
#include "Windows/Application.h"
#include "Graphics/EDXGui.h"
#include "Graphics/OpenGL.h"

#include "Base/Fluid.h"
#include "Smoke/Smoke.h"
#include "Water/Liquid.h"
#include "Util/VectorVisualizer.h"

using namespace EDX;
using namespace EDX::FluidSim;
using namespace EDX::GUI;

LiquidSolver<2> gFluid;
int giDim = 192;
const float gfTimeStep = 1.0f / 30.0f;
VectorFieldVisualizer2D gVisual;
bool gbSimulate = true;
int debugMode = 0;

void OnInit(Object* pSender, EventArgs args)
{
	OpenGL::InitializeOpenGLExtensions();

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

	gFluid.Initialize(Vector2i(giDim, giDim), Vector2::UNIT_SCALE);
	gVisual.Initialize(giDim, giDim, 14096);

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
	glPointSize(1.0f);
	ParticleVisualizer<2>::Render(gFluid.GetParticles(), 640.0f / float(giDim));
}

void OnResize(Object* pSender, ResizeEventArgs args)
{
	glViewport(0, 0, args.Width, args.Height);

	// Set opengl params
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, args.Width, 0, args.Height, -1, 1);

	glMatrixMode(GL_MODELVIEW);

	EDXGui::Resize(args.Width, args.Height);
}

void OnMouseEvent(Object* pSender, MouseEventArgs args)
{
	EDXGui::HandleMouseEvent(args);
}

void OnKeyboardEvent(Object* pSender, KeyboardEventArgs args)
{
	EDXGui::HandleKeyboardEvent(args);
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