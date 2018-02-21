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

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkAxesActor.h>
#include <vtkTransform.h>
#include <vtkSmartPointer.h>

#include "visu3d.h"

template<typename T>
using vSP = vtkSmartPointer<T>; // Shortcut

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
	vSP<vtkSphereSource> sphereSource = vSP<vtkSphereSource>::New();
	sphereSource->SetCenter(0.0, 0.0, 0.0);
	sphereSource->SetRadius(0.5);

	  //create a mapper
	vSP<vtkPolyDataMapper> sphereMapper = vSP<vtkPolyDataMapper>::New();
	sphereMapper->SetInputConnection(sphereSource->GetOutputPort());

	// create an actor
	vSP<vtkActor> sphereActor = vSP<vtkActor>::New();
	sphereActor->SetMapper(sphereMapper);

	// a renderer and render window
	vSP<vtkRenderer> renderer = vSP<vtkRenderer>::New();
	vSP<vtkRenderWindow> renderWindow = vSP<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);

	// an interactor
	vSP<vtkRenderWindowInteractor> renderWindowInteractor =
		vSP<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// add the actors to the scene
	renderer->AddActor(sphereActor);
	renderer->SetBackground(.1,.2,.3); // Background dark blue

	vSP<vtkTransform> transform = vSP<vtkTransform>::New();
	transform->Translate(1.0, 0.0, 0.0);

	vSP<vtkAxesActor> axes = vSP<vtkAxesActor>::New();

	// The axes are positioned with a user transform
	axes->SetUserTransform(transform);
	renderer->AddActor(axes);

	renderer->ResetCamera();
	renderWindow->Render();

	// begin mouse interaction
	renderWindowInteractor->Start();
}

