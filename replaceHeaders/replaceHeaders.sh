#!/bin/bash

find ./ -type f -name "*.C" -exec sed -i 1,22d {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    You should have received a copy of the GNU General Public License' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \\' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    for more details.' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \\' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    (at your option) any later version.' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    the Free Software Foundation, either version 3 of the License, or' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    under the terms of the GNU General Public License as published by' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    OpenFOAM is free software: you can redistribute it and/or modify it' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \\' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    This file is part of OpenFOAM.' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \License' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \-------------------------------------------------------------------------------' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    Copyright (C) 2019-2021 OpenCFD Ltd.' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    Copyright (C) 2011-2017 OpenFOAM Foundation' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \-------------------------------------------------------------------------------' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \     \\/     M anipulation  |' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \    \\  /    A nd           | www.openfoam.com' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \   \\    /   O peration     |' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \  =========                 |' {} +
find ./ -type f -name "*.C" -exec sed  -i '1i \/*---------------------------------------------------------------------------* \\' {} +


find ./ -type f -name "*.H" -exec sed -i 1,22d {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    You should have received a copy of the GNU General Public License' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \\' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    for more details.' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \\' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    (at your option) any later version.' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    the Free Software Foundation, either version 3 of the License, or' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    under the terms of the GNU General Public License as published by' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    OpenFOAM is free software: you can redistribute it and/or modify it' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \\' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    This file is part of OpenFOAM.' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \License' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \-------------------------------------------------------------------------------' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    Copyright (C) 2019-2021 OpenCFD Ltd.' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    Copyright (C) 2011-2017 OpenFOAM Foundation' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \-------------------------------------------------------------------------------' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \     \\/     M anipulation  |' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \    \\  /    A nd           | www.openfoam.com' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \   \\    /   O peration     |' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \  =========                 |' {} +
find ./ -type f -name "*.H" -exec sed  -i '1i \/*---------------------------------------------------------------------------* \\' {} +
