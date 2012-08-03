//reading the dictionary for the layer properties.
     scalar nLayers=readScalar(meshDict.lookup("nLayers"));
nLayers++;
//Calculating the length of the cell that's going to be refined (in weighted terms).
#    include "setRootCase.H"
#    include "createTime.H"
     runTime.functionObjects().off();
#    include "createPolyMesh.H"
     const word oldInstance = mesh.pointsInstance();
     word patchName
   (
         meshDict.lookup("patch")
     );
     bool overwrite = args.optionFound("overwrite");


    label patchID = mesh.boundaryMesh().findPatchID(patchName);

    if (patchID == -1)
    {
        FatalErrorIn(args.executable())
            << "Cannot find patch " << patchName << endl
            << "Valid patches are " << mesh.boundaryMesh().names()
            << exit(FatalError);
    }
    const polyPatch& pp = mesh.boundaryMesh()[patchID];


    // Cells cut

    labelHashSet cutCells(4*pp.size());

    const labelList& meshPoints = pp.meshPoints();

    forAll(meshPoints, pointI)
    {
        label meshPointI = meshPoints[pointI];

        const labelList& pCells = mesh.pointCells()[meshPointI];

        forAll(pCells, pCellI)
        {
            cutCells.insert(pCells[pCellI]);
        }
    }

    Info<< "Selected " << cutCells.size()
        << " cells connected to patch " << pp.name() << endl << endl;

    //
    // List of cells to refine
    //

    bool useSet = args.optionFound("useSet");

    if (useSet)
    {
        word setName(args.option("useSet"));

        Info<< "Subsetting cells to cut based on cellSet" << setName << endl
            << endl;

        cellSet cells(mesh, setName);

        Info<< "Read " << cells.size() << " cells from cellSet "
            << cells.instance()/cells.local()/cells.name()
            << endl << endl;

        for
        (
            cellSet::const_iterator iter = cells.begin();
            iter != cells.end();
            ++iter
        )
        {
            cutCells.erase(iter.key());
        }
        Info<< "Removed from cells to cut all the ones not in set " << setName
            << endl << endl;
    }

    // Mark all meshpoints on patch

    boolList vertOnPatch(mesh.nPoints(), false);

    forAll(meshPoints, pointI)
    {
        label meshPointI = meshPoints[pointI];

        vertOnPatch[meshPointI] = true;
    }


    // Mark cut edges.

    DynamicList<label> allCutEdges(pp.nEdges());

    DynamicList<scalar> allCutEdgeWeights(pp.nEdges());

    forAll(meshPoints, pointI)
    {
        label meshPointI = meshPoints[pointI];

        const labelList& pEdges = mesh.pointEdges()[meshPointI];

        forAll(pEdges, pEdgeI)
        {
            label edgeI = pEdges[pEdgeI];

            const edge& e = mesh.edges()[edgeI];

            label otherPointI = e.otherVertex(meshPointI);

            if (!vertOnPatch[otherPointI])
            {
                allCutEdges.append(edgeI);

                if (e.start() == meshPointI)
                {
                    allCutEdgeWeights.append(weight);
                }
                else
                {
                    allCutEdgeWeights.append(1 - weight);
                }
            }
        }
    }

    allCutEdges.shrink();
    allCutEdgeWeights.shrink();

    Info<< "Cutting:" << endl
        << "    cells:" << cutCells.size() << endl
        << "    edges:" << allCutEdges.size() << endl
        << endl;

    // Transfer DynamicLists to straight ones.
    scalarField cutEdgeWeights;
    cutEdgeWeights.transfer(allCutEdgeWeights);
    allCutEdgeWeights.clear();


    // Gets cuts across cells from cuts through edges.
    cellCuts cuts
    (
        mesh,
        cutCells.toc(),     // cells candidate for cutting
        labelList(0),       // cut vertices
        allCutEdges,        // cut edges
        cutEdgeWeights      // weight on cut edges
    );

    directTopoChange meshMod(mesh);

    // Cutting engine
    meshCutter cutter(mesh);

    // Insert mesh refinement into directTopoChange.
    cutter.setRefinement(cuts, meshMod);

    // Do all changes
    Info<< "Morphing ..." << endl;

    if (!overwrite)
    {
        runTime++;
    }

    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    if (morphMap().hasMotionPoints())
    {
        mesh.movePoints(morphMap().preMotionPoints());
    }

    // Update stored labels on meshCutter.
    cutter.updateMesh(morphMap());

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    // Write resulting mesh
    Info << "Writing refined morphMesh to time " << runTime.timeName() << endl;

    mesh.write();

    Info << "End\n" << endl;




// ************************************************************************* //
