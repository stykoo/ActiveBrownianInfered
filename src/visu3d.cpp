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

#include <vector>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
//#include <vtkAxesActor.h>
#include <vtkTransform.h>
#include <vtkSmartPointer.h>
#include <vtkCommand.h>

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

			// Update the positions
			for (size_t i = 0 ; i < sphereActors->size() ; ++i) {
				(*sphereActors)[i]->SetPosition(state->getPosX(i),
                                                state->getPosY(i),
											    state->getPosZ(i));
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

	for (long i = 0 ; i < n_parts ; ++i) {
		sphereSources.push_back(vSP<vtkSphereSource>::New());
		//sphereSources[i]->SetCenter(state->getPosX(i), state->getPosY(i),
				                    //state->getPosZ(i));
		sphereSources[i]->SetRadius(0.5);

		//create a mapper
		sphereMappers.push_back(vSP<vtkPolyDataMapper>::New());
		sphereMappers[i]->SetInputConnection(sphereSources[i]->GetOutputPort());

		// create an actor
		sphereActors.push_back(vSP<vtkActor>::New());
		sphereActors[i]->SetMapper(sphereMappers[i]);
		sphereActors[i]->SetPosition(state->getPosX(i), state->getPosY(i),
				                     state->getPosZ(i));

	}

	// a renderer and render window
	vSP<vtkRenderer> renderer = vSP<vtkRenderer>::New();
	vSP<vtkRenderWindow> renderWindow = vSP<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);

	// an interactor
	vSP<vtkRenderWindowInteractor> renderWindowInteractor =
		vSP<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// add the actors to the scene
	for (long i = 0 ; i < n_parts ; ++i) {
		renderer->AddActor(sphereActors[i]);
	}
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
