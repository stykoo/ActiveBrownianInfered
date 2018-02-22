/*
Copyright (C) Sorbonne Université (2018)
Contributor: Alexis Poncet <aponcet@lptmc.jussieu.fr>

This file is part of ActiveBrownian.

ActiveBrownian is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ActiveBrownian is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ActiveBrownian.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file visu3d.cpp
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Visualization of the system in 3d.
 *
 * The system is visualized using the VTK library.
*/

#include <iostream>
#include <vector>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkCubeSource.h>
#include <vtkAxesActor.h>
#include <vtkTransform.h>
#include <vtkSmartPointer.h>
#include <vtkCommand.h>
#include <vtkProperty.h>

#include "visu3d.h"

template<typename T>
using vSP = vtkSmartPointer<T>; // Shortcut

struct Visu3dTimer : public vtkCommand {
	public:
		vtkTypeMacro(Visu3dTimer, vtkCommand);
		static Visu3dTimer *New() {
    		return new Visu3dTimer();
  		}

		void Execute(vtkObject *caller, unsigned long vtkNotUsed(eventId),
					   void *vtkNotUsed(callData)) {
			vtkRenderWindowInteractor *iren =
				static_cast<vtkRenderWindowInteractor*>(caller);

			// Update the positions and the colors
			double r, g, b;
			for (size_t i = 0 ; i < sphereActors->size() ; ++i) {
				(*sphereActors)[i]->SetPosition(state->getPosX(i),
                                                state->getPosY(i),
											    state->getPosZ(i));
				colorFromAngles(r, g, b, state->getOrient(i)->getTheta(),
								state->getOrient(i)->getPhi()); 
				(*sphereActors)[i]->GetProperty()->SetColor(r, g, b);
			}

			iren->Render();
		}

	// Public attributes (this is dirty...)
	std::vector< vSP<vtkActor> > *sphereActors;
	const State3d *state;
};

/*!
 * \brief Constructor for visualization
 *
 * \param state Pointer to the state of the system
 * \param len Length of the box
 * \param n_parts Number of particles
 */
Visu3d::Visu3d(const State3d *state, const double len, const long n_parts) :
	state(state), len(len), n_parts(n_parts) {
}

/*!
 * \brief Thread for visualization.
 *
 * Open a window, draw the particles and update their
 * positions while the simulation is runing.
 */
void Visu3d::run() {
	std::vector< vSP<vtkSphereSource> > sphereSources;
	std::vector< vSP<vtkPolyDataMapper> > sphereMappers;
	std::vector< vSP<vtkActor> > sphereActors;
	double r, g, b;

	for (long i = 0 ; i < n_parts ; ++i) {
		sphereSources.push_back(vSP<vtkSphereSource>::New());
		sphereSources[i]->SetRadius(0.5);
		sphereSources[i]->SetThetaResolution(sphere_res);
		sphereSources[i]->SetPhiResolution(sphere_res);

		//create a mapper
		sphereMappers.push_back(vSP<vtkPolyDataMapper>::New());
		sphereMappers[i]->SetInputConnection(sphereSources[i]->GetOutputPort());

		// create an actor
		sphereActors.push_back(vSP<vtkActor>::New());
		sphereActors[i]->SetMapper(sphereMappers[i]);
		sphereActors[i]->SetPosition(state->getPosX(i), state->getPosY(i),
				                     state->getPosZ(i));
		colorFromAngles(r, g, b, state->getOrient(i)->getTheta(),
				        state->getOrient(i)->getPhi()); 
		sphereActors[i]->GetProperty()->SetColor(r, g, b);
		sphereActors[i]->GetProperty()->SetOpacity(sphere_opa);

	}

	// Box
	vSP<vtkCubeSource> cubeSource = vSP<vtkCubeSource>::New();
	cubeSource->SetXLength(len);
	cubeSource->SetYLength(len);
	cubeSource->SetZLength(len);
	cubeSource->SetCenter(len / 2., len / 2., len / 2.);
	cubeSource->Update();
	vSP<vtkPolyDataMapper> cubeMapper = vSP<vtkPolyDataMapper>::New();
	cubeMapper->SetInputConnection(cubeSource->GetOutputPort());
	vSP<vtkActor> cubeActor = vSP<vtkActor>::New();
	cubeActor->GetProperty()->SetRepresentationToWireframe();
	cubeActor->SetMapper(cubeMapper);

	// Axes
	vSP<vtkAxesActor> axes = vSP<vtkAxesActor>::New();
	axes->SetXAxisLabelText("");
	axes->SetYAxisLabelText("");
	axes->SetZAxisLabelText("");

	// Renderer and render window
	vSP<vtkRenderer> renderer = vSP<vtkRenderer>::New();
	vSP<vtkRenderWindow> renderWindow = vSP<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);

	// Interactor
	vSP<vtkRenderWindowInteractor> renderWindowInteractor =
		vSP<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Add the actors to the scene
	for (long i = 0 ; i < n_parts ; ++i) {
		renderer->AddActor(sphereActors[i]);
	}
	renderer->AddActor(cubeActor);
	renderer->AddActor(axes);
	renderer->SetBackground(.1,.2,.3); // Background dark blue

	// Take care of animation
	renderWindowInteractor->Initialize();
	renderWindowInteractor->CreateRepeatingTimer(delay);
	vtkSmartPointer<Visu3dTimer> timerCallback =
		vtkSmartPointer<Visu3dTimer>::New();
	timerCallback->sphereActors = &sphereActors;
	timerCallback->state = state;
	renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, timerCallback);

	// Render and interact
	renderWindow->Render();
	renderWindowInteractor->Start();
}

/*!
 * \brief Associate a color to a set of angles (theta, phi)
 *
 * Use HSL representation with H = phi in degrees, S = 1, V = theta / pi  
 * 
 * \param rgb Array of RGB color to be modfied
 * \param theta Angle between 0 and pi
 * \param phi Angle between 0 and 2pi
*/
void colorFromAngles(double &r, double &g, double &b,
		             const double theta, const double phi) {
	/*double x = std::sin(theta) * std::cos(phi);
	double y = std::sin(theta) * std::sin(phi);
	double z = std::cos(theta);
	r = 0.5 * (1 + x);
	g = 0.5 * (1 + y);
	b = 0.5 * (1 + z); */
	int hue = (int) (phi * 180. / M_PI) + 180;
	double light = 0.5 + 0.9 * (theta / M_PI - 0.5);
	double sat_hsl = 1.0;

	// HSL to HSV
	double a = sat_hsl * (light<0.5 ? light : 1-light);
	double val = light + a;
  	double sat = light>0 ? 2*a/val : 0;

	int h = hue/60;
	double f = double(hue)/60-h;
	double p = val*(1.f-sat);
	double q = val*(1.f-sat*f);
	double t = val*(1.f-sat*(1-f));

	switch(h)
	{
		default:
		case 0:
		case 6:
			r = val; g = t; b = p;
			break;
		case 1:
			r = q; g = val; b = p;
			break;
		case 2:
			r = p; g = val; b = t;
			break;
		case 3:
		 	r = p; g = q; b = val;
			break;
		case 4:
			r = t; g = p; b = val;
			break;
		case 5:
			r = val; g = p; b = q;
			break;
	}
}
