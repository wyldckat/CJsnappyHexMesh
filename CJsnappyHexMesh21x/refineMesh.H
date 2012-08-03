#   include "createPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    printEdgeStats(mesh);

pointMesh pMesh(mesh);
pointVectorField pointMotionU
(
     IOobject
  (
          "pointMotionU",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
  ),
      pMesh
);


    //
    // Read/construct control dictionary
    //

    bool readDict = args.optionFound("dict");
    bool overwrite = args.optionFound("overwrite");

    // List of cells to refine
    labelList refCells;

    // Dictionary to control refinement
    dictionary refineDict;

    if (readDict)
    {
        Info<< "Refining according to refineMeshDict" << nl << endl;

        refineDict =
            IOdictionary
            (
                IOobject
                (
                    "refineMeshDict",
                    runTime.system(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

        word setName(refineDict.lookup("set"));

        cellSet cells(mesh, setName);

        Pout<< "Read " << cells.size() << " cells from cellSet "
            << cells.instance()/cells.local()/cells.name()
            << endl << endl;

        refCells = cells.toc();
    }
    else
    {
        Info<< "Refining all cells" << nl << endl;

        // Select all cells
        refCells.setSize(mesh.nCells());

        forAll(mesh.cells(), cellI)
        {
            refCells[cellI] = cellI;
        }


        // Set refinement directions based on 2D/3D
        label axisIndex = twoDNess(mesh);

        if (axisIndex == -1)
        {
            Info<< "3D case; refining all directions" << nl << endl;

            wordList directions(3);
            directions[0] = "tan1";
            directions[1] = "tan2";
            directions[2] = "normal";
            refineDict.add("directions", directions);

            // Use hex cutter
            refineDict.add("useHexTopology", "true");
        }
        else
        {
            wordList directions(2);

            if (axisIndex == 0)
            {
                Info<< "2D case; refining in directions y,z\n" << endl;
                directions[0] = "tan2";
                directions[1] = "normal";
            }
            else if (axisIndex == 1)
            {
                Info<< "2D case; refining in directions x,z\n" << endl;
                directions[0] = "tan1";
                directions[1] = "normal";
            }
            else
            {
                Info<< "2D case; refining in directions x,y\n" << endl;
                directions[0] = "tan1";
                directions[1] = "tan2";
            }

            refineDict.add("directions", directions);

            // Use standard cutter
            refineDict.add("useHexTopology", "false");
        }

        refineDict.add("coordinateSystem", "global");

        dictionary coeffsDict;
        coeffsDict.add("tan1", vector(1, 0, 0));
        coeffsDict.add("tan2", vector(0, 1, 0));
        refineDict.add("globalCoeffs", coeffsDict);

        refineDict.add("geometricCut", "false");
        refineDict.add("writeMesh", "false");
    }


    string oldTimeName(runTime.timeName());

    if (!overwrite)
    {
        runTime++;
    }


    // Multi-directional refinement (does multiple iterations)
    multiDirRefinement multiRef(mesh, refCells, refineDict);


    // Write resulting mesh
    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    mesh.write();


    // Get list of cell splits.
    // (is for every cell in old mesh the cells they have been split into)
    const labelListList& oldToNew = multiRef.addedCells();


    // Create cellSet with added cells for easy inspection
    cellSet newCells(mesh, "refinedCells", refCells.size());

    forAll(oldToNew, oldCellI)
    {
        const labelList& added = oldToNew[oldCellI];

        forAll(added, i)
        {
            newCells.insert(added[i]);
        }
    }

    Pout<< "Writing refined cells (" << newCells.size() << ") to cellSet "
        << newCells.instance()/newCells.local()/newCells.name()
        << endl << endl;

    newCells.write();




    //
    // Invert cell split to construct map from new to old
    //

    labelIOList newToOld
    (
        IOobject
        (
            "cellMap",
            runTime.timeName(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.nCells()
    );
    newToOld.note() =
        "From cells in mesh at "
      + runTime.timeName()
      + " to cells in mesh at "
      + oldTimeName;


    forAll(oldToNew, oldCellI)
    {
        const labelList& added = oldToNew[oldCellI];

        if (added.size())
        {
            forAll(added, i)
            {
                newToOld[added[i]] = oldCellI;
            }
        }
        else
        {
            // Unrefined cell
            newToOld[oldCellI] = oldCellI;
        }
    }

    Info<< "Writing map from new to old cell to "
        << newToOld.objectPath() << nl << endl;

    newToOld.write();


    // Some statistics.

    printEdgeStats(mesh);

    Info<< "End\n" << endl;





// ************************************************************************* //
