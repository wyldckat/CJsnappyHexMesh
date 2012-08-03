#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"

    for(int i=0;i<nDynIter;i++)
    {
        if (!overwrite)
        {
            const_cast<Time&>(mesh.time())++;
        }
    	Info<< "Time = " << runTime.timeName() << endl;

        if (isDir(runTime.path()/"VTK"))
        {
            Info << "Clear VTK directory" << endl;
            rmDir(runTime.path()/"VTK");
        }

        mesh.update();

//#       include "checkVolContinuity.H"

        mesh.checkMesh(true);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;




// ************************************************************************* //
