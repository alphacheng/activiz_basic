using System;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using System.Drawing;
using System.Diagnostics;

using Kitware.VTK;
using static System.Net.Mime.MediaTypeNames;
using System.Reflection;


namespace ActivizLearn
{
    public partial class Form1 : Form
    {
        private RenderWindowControl myRenderWindowControl;
        public Form1()
        {
            InitializeComponent();
            myRenderWindowControl = new RenderWindowControl();
            myRenderWindowControl.SetBounds(0, 0, panel1.Width, panel1.Height);
            myRenderWindowControl.Dock = DockStyle.Fill;
            panel1.Controls.Add(this.myRenderWindowControl);

            myRenderWindowControl.Load += RenderWindowControl1_Load;
        }

        private void RenderWindowControl1_Load(object sender, System.EventArgs e)
        {
            try
            {
                //DrawPoint();
                //DrawTriangle();
                //DrawCylinder();
                //DrawSphere(0.5);
                //ReadStlFileDraw();
                //DrawAssembly();
                //DrawTexturePlane();
                //DrawRainBow();
                //DrawVisQuad();
                //DrawBuildUGrid();

                //PolygonalSurfaceContourLineInterpolator();

                //vtkPolyDataConnectivityFilter_LargestRegion();

                //WeightedTransformFilter();

                //Decimation();

                //Subdivision();

                //ExtractEdges();

                //ColoredElevationMap();

                //MarchingCubes();

                //ExtractLargestIsoSurface("test.vtk", 50, true);
                //test.vtk: http://web.kaist.ac.kr/~hjy/test.vtk

                //VTKSurfaceReconstruction();

                polyDataNormals();








            }
            catch(Exception ex)
            {
                MessageBox.Show(ex.Message, "Exception", MessageBoxButtons.OK);
            }
        }


        private void polyDataNormals()
        {
            vtkPolyDataReader plyReader =
                  vtkPolyDataReader.New();
            plyReader.SetFileName("../../data/fran_cut.vtk");
            plyReader.Update();

            vtkPolyDataNormals normFilter =
                  vtkPolyDataNormals.New();
            normFilter.SetInput(plyReader.GetOutput());
            normFilter.SetComputePointNormals(1);//开启点法向量计算
            normFilter.SetComputeCellNormals(0); //关闭单元法向量计算
            normFilter.SetAutoOrientNormals(1);
            normFilter.SetSplitting(0);
            normFilter.Update();

            vtkMaskPoints mask =
                  vtkMaskPoints.New();
            mask.SetInput(normFilter.GetOutput());
            mask.SetMaximumNumberOfPoints(300);
            mask.RandomModeOn();
            mask.Update();

            vtkArrowSource arrow =
                  vtkArrowSource.New();
            arrow.Update(); //一定要更新 否则数据没有添加进来，程序会报错

            vtkGlyph3D glyph =
                  vtkGlyph3D.New();
            glyph.SetInput(mask.GetOutput());
            glyph.SetSource(arrow.GetOutput());//每一点用箭头代替
            glyph.SetVectorModeToUseNormal();//设置向量显示模式和法向量一致
            glyph.SetScaleFactor(0.01); //设置伸缩比例
            glyph.Update();

            vtkPolyDataMapper mapper =
                  vtkPolyDataMapper.New();
            mapper.SetInput(plyReader.GetOutput());
            vtkPolyDataMapper normMapper =
                  vtkPolyDataMapper.New();
            normMapper.SetInput(normFilter.GetOutput());
            vtkPolyDataMapper glyphMapper =
                  vtkPolyDataMapper.New();
            glyphMapper.SetInput(glyph.GetOutput());

            vtkActor actor =
                  vtkActor.New();
            actor.SetMapper(mapper);
            vtkActor normActor =
                  vtkActor.New();
            normActor.SetMapper(normMapper);
            vtkActor glyphActor =
                  vtkActor.New();
            glyphActor.SetMapper(glyphMapper);
            glyphActor.GetProperty().SetColor(1, 0, 0);

            double[] origView = new double[] { 0, 0, 0.33, 1 };
            double[] normView = new double[] { 0.33, 0, 0.66, 1 };
            double[] glyphView = new double[] { 0.66, 0, 1, 1 };
            vtkRenderer origRender =  vtkRenderer.New();
            origRender.SetViewport(origView[0], origView[1], origView[2], origView[3]);
            origRender.AddActor(actor);
            origRender.SetBackground(1, 0, 0);
            vtkRenderer normRender =  vtkRenderer.New();
            normRender.SetViewport(normView[0], normView[1], normView[2], normView[3]);
            normRender.AddActor(normActor);
            normRender.SetBackground(0, 1, 0);
            vtkRenderer glyphRender = vtkRenderer.New();
            glyphRender.SetViewport(glyphView[0], glyphView[1], glyphView[2], glyphView[3]);
            glyphRender.AddActor(glyphActor);
            glyphRender.AddActor(normActor);
            glyphRender.SetBackground(0, 0, 1);

            // get a reference to the renderwindow of our renderWindowControl1
            vtkRenderWindow renderWindow = myRenderWindowControl.RenderWindow;
            renderWindow.AddRenderer(origRender);
            renderWindow.AddRenderer(normRender);
            renderWindow.AddRenderer(glyphRender);
            renderWindow.SetWindowName("Calculating Point Norm & Cell Norm");
            renderWindow.SetSize(960, 320);
            renderWindow.Render();


        }
        private void VTKSurfaceReconstruction()
        {
            vtkPolyDataReader reader = vtkPolyDataReader.New();

            Debug.WriteLine(Environment.CurrentDirectory);

            reader.SetFileName("../../data/fran_cut.vtk");
            reader.Update();

        vtkPolyData points =vtkPolyData.New();
        points.SetPoints(reader.GetOutput().GetPoints()); //获得网格模型中的几何数据：点集

        vtkSurfaceReconstructionFilter surf = vtkSurfaceReconstructionFilter.New();
        surf.SetInput(points);
        surf.SetNeighborhoodSize(20);
        surf.SetSampleSpacing(0.005);
        surf.Update();

        vtkContourFilter contour = vtkContourFilter.New();
        contour.SetInputConnection(surf.GetOutputPort());
        contour.SetValue(0, 0.0);
        contour.Update();
        //
        vtkVertexGlyphFilter vertexGlyphFilter = vtkVertexGlyphFilter.New();
        vertexGlyphFilter.AddInput(points);
        vertexGlyphFilter.Update();
        vtkPolyDataMapper pointMapper = vtkPolyDataMapper.New();
        pointMapper.SetInput(vertexGlyphFilter.GetOutput());
        pointMapper.ScalarVisibilityOff();

        vtkActor pointActor = vtkActor.New();
        pointActor.SetMapper(pointMapper);
        pointActor.GetProperty().SetColor(1, 0, 0);
        pointActor.GetProperty().SetPointSize(4);

        vtkPolyDataMapper contourMapper = vtkPolyDataMapper.New();
        contourMapper.SetInput(contour.GetOutput());

        for (int i = 0; i<3; i++)
            {
                Debug.WriteLine((contour.GetOutput().GetBounds()[i]- contour.GetOutput().GetBounds()[i+1]).ToString() +"vs" 
                + (vertexGlyphFilter.GetOutput().GetBounds()[i]- vertexGlyphFilter.GetOutput().GetBounds()[i+1]).ToString());
            }


        vtkActor contourActor = vtkActor.New();
        contourActor.SetMapper(contourMapper);
        ///

        double[] pointView = new double[] { 0.0, 0.0, 0.5, 1.0 };
        double[] contourView = new double[] { 0.5, 0.0, 1.0, 1.0 };

        vtkRenderer pointRender = vtkRenderer.New();
        pointRender.AddActor(pointActor);
        pointRender.SetViewport(pointView[0], pointView[1], pointView[2], pointView[3]);
        pointRender.SetBackground(1, 1, 1);

        vtkRenderer contourRender = vtkRenderer.New();
        contourRender.AddActor(contourActor);
        //contourRender.AddActor(pointActor);
        contourRender.SetViewport(contourView[0], contourView[1], contourView[2], contourView[3]);
        contourRender.SetBackground(0, 1, 0);

            vtkCamera camera = vtkCamera.New();
            camera.SetPosition(0, -1, 0);
            camera.SetFocalPoint(0, 0, 0);
            camera.SetViewUp(0, 0, 1);
            camera.Azimuth(30);
            camera.Elevation(30);
                 //pointRender.ResetCamera();
        pointRender.SetActiveCamera(camera);
        contourRender.SetActiveCamera(camera);

            // get a reference to the renderwindow of our renderWindowControl1
            vtkRenderWindow renderWindow = myRenderWindowControl.RenderWindow;
            renderWindow.AddRenderer(pointRender);
            renderWindow.AddRenderer(contourRender);
            renderWindow.SetSize(640, 320);
            renderWindow.SetWindowName("3D Surface Reconstruction ");
            renderWindow.Render();


        }

        private void ExtractLargestIsoSurface(string filePath, double threshold, bool extractLargest)
        {
            // Load data
            vtkStructuredPointsReader reader = vtkStructuredPointsReader.New();
            reader.SetFileName(filePath);

            // Create a 3D model using marching cubes
            vtkMarchingCubes mc = vtkMarchingCubes.New();
            mc.SetInputConnection(reader.GetOutputPort());
            mc.ComputeNormalsOn();
            mc.ComputeGradientsOn();
            mc.SetValue(0, threshold);  // second value acts as threshold

            // To remain largest region
            vtkPolyDataConnectivityFilter confilter = vtkPolyDataConnectivityFilter.New();
            confilter.SetInputConnection(mc.GetOutputPort());
            confilter.SetExtractionModeToLargestRegion();

            // Create a mapper
            vtkPolyDataMapper mapper = vtkPolyDataMapper.New();
            if (extractLargest)
            {
                mapper.SetInputConnection(confilter.GetOutputPort());
            }
            else
            {
                mapper.SetInputConnection(mc.GetOutputPort());
            }

            mapper.ScalarVisibilityOff();    // utilize actor's property I set

            // Visualize
            vtkActor actor = vtkActor.New();
            actor.GetProperty().SetColor(1, 1, 1);
            actor.SetMapper(mapper);

            // get a reference to the renderwindow of our renderWindowControl1
            vtkRenderWindow renderWindow = myRenderWindowControl.RenderWindow;
            // renderer
            vtkRenderer renderer = renderWindow.GetRenderers().GetFirstRenderer();
            // add our actor to the renderer
            renderer.AddActor(actor);
        }

        private void MarchingCubes()
        {
            vtkSphereSource sphereSource = vtkSphereSource.New();
            sphereSource.SetPhiResolution(20);
            sphereSource.SetThetaResolution(20);
            sphereSource.Update();

            double[] bounds = sphereSource.GetOutput().GetBounds();
            for (int i = 0; i<6; i += 2)
            {
                double range = bounds[i + 1] - bounds[i];
                bounds[i] = bounds[i] - .1 * range;
                bounds[i + 1] = bounds[i + 1] + .1 * range;
            }
            vtkVoxelModeller voxelModeller = vtkVoxelModeller.New();
            voxelModeller.SetSampleDimensions(50, 50, 50);
            voxelModeller.SetModelBounds(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
            voxelModeller.SetScalarTypeToFloat();
            voxelModeller.SetMaximumDistance(.1);

#if VTK_MAJOR_VERSION_5
         voxelModeller.SetInputConnection(sphereSource.GetOutputPort());
#else
            voxelModeller.SetInput(sphereSource.GetOutput());
#endif
            vtkMarchingCubes surface = vtkMarchingCubes.New();

#if VTK_MAJOR_VERSION_5
         surface.SetInputConnection(voxelModeller.GetOutputPort());
#else
            surface.SetInput(voxelModeller.GetOutput());
#endif
            surface.ComputeNormalsOn();
            surface.SetValue(0, 0.5);
            vtkPolyDataMapper mapper = vtkPolyDataMapper.New();
#if VTK_MAJOR_VERSION_5
         mapper.SetInputConnection(surface.GetOutputPort());
#else
            mapper.SetInput(surface.GetOutput());
#endif
            vtkActor actor = vtkActor.New();
            actor.SetMapper(mapper);

            // get a reference to the renderwindow of our renderWindowControl1
            vtkRenderWindow renderWindow = myRenderWindowControl.RenderWindow;
            // renderer
            vtkRenderer renderer = renderWindow.GetRenderers().GetFirstRenderer();
            // set background color
            renderer.SetBackground(.2, .3, .4);
            // add our actor to the renderer
            renderer.AddActor(actor);
        }


        private void ColoredElevationMap()
        {
            // Create a grid of points (height/terrian map)
            vtkPoints points = vtkPoints.New();

            uint GridSize = 20;
            double xx, yy, zz;
            for (uint x = 0; x<GridSize; x++)
            {
                for (uint y = 0; y<GridSize; y++)
                {
                    xx = x + vtkMath.Random(-.2, .2);
                    yy = y + vtkMath.Random(-.2, .2);
                    zz = vtkMath.Random(-.5, .5);
                    points.InsertNextPoint(xx, yy, zz);
                }
            }

            // Add the grid points to a polydata object
            vtkPolyData inputPolyData = vtkPolyData.New();
            inputPolyData.SetPoints(points);

            // Triangulate the grid points
            vtkDelaunay2D delaunay = vtkDelaunay2D.New();
#if VTK_MAJOR_VERSION_5
         delaunay.SetInput(inputPolyData);
#else
            delaunay.SetInput(inputPolyData);
#endif
            delaunay.Update();
            vtkPolyData outputPolyData = delaunay.GetOutput();
            vtkPolyData outputPolyData2 = delaunay.GetOutput();

            double[] bounds = outputPolyData.GetBounds();

            // Find min and max z
            double minz = bounds[4];
            double maxz = bounds[5];

            Debug.WriteLine("minz: " + minz);
            Debug.WriteLine("maxz: " + maxz);

            // Create the color map
            vtkLookupTable colorLookupTable = vtkLookupTable.New();
            colorLookupTable.SetTableRange(minz, maxz);
            colorLookupTable.SetNumberOfTableValues(9);
            colorLookupTable.Build();

            // Create the color map
            vtkLookupTable colorLookupTable2 = vtkLookupTable.New();
            colorLookupTable2.SetTableRange(minz, maxz);
            //colorLookupTable.SetNumberOfTableValues(9);
            colorLookupTable2.Build();

            // Generate the colors for each point based on the color map
            vtkUnsignedCharArray colors = vtkUnsignedCharArray.New();
            colors.SetNumberOfComponents(3);
            colors.SetName("Colors");

            // Generate the colors for each point based on the color map
            vtkUnsignedCharArray colors2 = vtkUnsignedCharArray.New();
            colors2.SetNumberOfComponents(3);
            colors2.SetName("Colors");

            Debug.WriteLine("There are " + outputPolyData.GetNumberOfPoints()
                      + " points.");


#if UNSAFE // fastest way to fill color array
         colors.SetNumberOfTuples(outputPolyData.GetNumberOfPoints());
         unsafe {
            byte* pColor = (byte*)colors.GetPointer(0).ToPointer();

            for(int i = 0; i  outputPolyData.GetNumberOfPoints(); i++) {
               double[] p = outputPolyData.GetPoint(i);

               double[] dcolor = colorLookupTable.GetColor(p[2]);
               Debug.WriteLine("dcolor: "
                         + dcolor[0] + " "
                         + dcolor[1] + " "
                         + dcolor[2]);

               byte[] color = new byte[3];
               for(uint j = 0; j  3; j++) {
                  color[j] = (byte)( 255 * dcolor[j] );
               }
               Debug.WriteLine("color: "
                         + color[0] + " "
                         + color[1] + " "
                         + color[2]);

               *( pColor + 3 * i ) = color[0];
               *( pColor + 3 * i + 1 ) = color[1];
               *( pColor + 3 * i + 2 ) = color[2];
            }
         }
#else
            for (int i = 0; i<outputPolyData.GetNumberOfPoints(); i++)
                 
            {
                double[] p = outputPolyData.GetPoint(i);

                double[] dcolor = colorLookupTable.GetColor(p[2]);
                double[] dcolor2 = colorLookupTable2.GetColor(p[2]);
                Debug.WriteLine("dcolor: "
                          + dcolor[0] + " "
                          + dcolor[1] + " "
                          + dcolor[2]);

                byte[] color = new byte[3];
                byte[] color2 = new byte[3];
                for (uint j = 0; j<3; j++)
                {
                    color[j] = (byte)(255 * dcolor[j]);
                    color2[j] = (byte)(255 * dcolor2[j]);
                }
                Debug.WriteLine("color: "
                          + color[0] + " "
                          + color[1] + " "
                          + color[2]);
                colors.InsertNextTuple3(color[0], color[1], color[2]);
                colors2.InsertNextTuple3(color2[0], color2[1], color2[2]);
                //IntPtr pColor = Marshal.AllocHGlobal(Marshal.SizeOf(typeof(byte)) * 3);
                //Marshal.Copy(color, 0, pColor, 3);
                //colors.InsertNextTupleValue(pColor);
                //Marshal.FreeHGlobal(pColor);
            }
#endif

            outputPolyData.GetPointData().SetScalars(colors);
            outputPolyData2.GetPointData().SetScalars(colors2);

            // Create a mapper and actor
            vtkPolyDataMapper mapper = vtkPolyDataMapper.New();
            vtkPolyDataMapper mapper2 = vtkPolyDataMapper.New();
#if VTK_MAJOR_VERSION_5
         mapper.SetInputConnection(outputPolyData.GetProducerPort());
#else
            mapper.SetInput(outputPolyData);
            mapper2.SetInput(outputPolyData2);
#endif

            vtkActor actor = vtkActor.New();
            actor.SetMapper(mapper);

            vtkActor actor2 = vtkActor.New();
            actor2.SetMapper(mapper2);

            // get a reference to the renderwindow of our renderWindowControl1
            vtkRenderWindow renderWindow = myRenderWindowControl.RenderWindow;



            // Define viewport ranges
            // (xmin, ymin, xmax, ymax)
            double[] leftViewport = new double[] { 0.0, 0.0, 0.5, 1.0 };
            double[] rightViewport = new double[] { 0.5, 0.0, 1.0, 1.0 };

            // Setup both renderers
            vtkRenderer leftRenderer = vtkRenderer.New();
            renderWindow.AddRenderer(leftRenderer);
            leftRenderer.SetViewport(leftViewport[0], leftViewport[1], leftViewport[2], leftViewport[3]);
            leftRenderer.SetBackground(.6, .5, .4);

            vtkRenderer rightRenderer = vtkRenderer.New();
            renderWindow.AddRenderer(rightRenderer);
            rightRenderer.SetViewport(rightViewport[0], rightViewport[1], rightViewport[2], rightViewport[3]);
            rightRenderer.SetBackground(.4, .5, .6);

            // Add the sphere to the left and the cube to the right
            leftRenderer.AddActor(actor);
            rightRenderer.AddActor(actor2);
            leftRenderer.ResetCamera();
            rightRenderer.ResetCamera();
            renderWindow.Render();
        }

        private void ExtractEdges()
        {
            vtkSphereSource sphereSource = vtkSphereSource.New();
            sphereSource.Update();

            Debug.WriteLine("Sphere" + Environment.NewLine + "----------");
            Debug.WriteLine("There are " + sphereSource.GetOutput().GetNumberOfCells() + " cells.");
            Debug.WriteLine("There are " + sphereSource.GetOutput().GetNumberOfPoints() + " points.");

            vtkExtractEdges extractEdges = vtkExtractEdges.New();
#if VTK_MAJOR_VERSION_5
         extractEdges.SetInputConnection(sphereSource.GetOutputPort());
#else
            extractEdges.SetInput(sphereSource.GetOutput());
#endif
            extractEdges.Update();

            vtkCellArray lines = extractEdges.GetOutput().GetLines();
            vtkPoints points = extractEdges.GetOutput().GetPoints();

            Debug.WriteLine(Environment.NewLine + "Edges" + Environment.NewLine + "----------");
            Debug.WriteLine("There are " + lines.GetNumberOfCells() + " cells.");
            Debug.WriteLine("There are " + points.GetNumberOfPoints() + " points.");

            // Traverse all of the edges
            for (int i = 0; i<extractEdges.GetOutput().GetNumberOfCells(); i++)
            {
                //Debug.WriteLine("Type: " + extractEdges.GetOutput().GetCell(i).GetClassName() );
                vtkLine line = vtkLine.SafeDownCast(extractEdges.GetOutput().GetCell(i));
                Debug.WriteLine("Line " + i + " : " + line);
            }

            // Visualize the edges

            // Create a mapper and actor
            vtkPolyDataMapper mapper = vtkPolyDataMapper.New();
#if VTK_MAJOR_VERSION_5
         mapper.SetInputConnection(extractEdges.GetOutputPort());
#else
            mapper.SetInput(extractEdges.GetOutput());
#endif
            vtkActor actor = vtkActor.New();
            actor.SetMapper(mapper);
            // get a reference to the renderwindow of our renderWindowControl1
            vtkRenderWindow renderWindow = myRenderWindowControl.RenderWindow;
            // renderer
            vtkRenderer renderer = renderWindow.GetRenderers().GetFirstRenderer();
            // set background color
            renderer.SetBackground(1, 1, 1);
            // add our actor to the renderer
            renderer.AddActor(actor);
        }


        private void Subdivision()
        {
            vtkSphereSource sphereSource = vtkSphereSource.New();
            sphereSource.Update();

            vtkPolyData input = vtkPolyData.New();
            input.ShallowCopy(sphereSource.GetOutput());

            Debug.WriteLine("Before decimation" + Environment.NewLine + "------------");
            Debug.WriteLine("There are " + input.GetNumberOfPoints() + " points.");
            Debug.WriteLine("There are " + input.GetNumberOfPolys() + " polygons.");

            vtkLoopSubdivisionFilter loop = vtkLoopSubdivisionFilter.New();
#if VTK_MAJOR_VERSION_5
         decimate.SetInputConnection(input.GetProducerPort());
#else
            loop.SetInput(input);
#endif
            //decimate.SetTargetReduction(.99); //99% reduction (if there was 100 triangles, now there will be 1)
            loop.SetNumberOfSubdivisions(1); //10% reduction (if there was 100 triangles, now there will be 50)
            loop.Update();

            vtkPolyData decimated = vtkPolyData.New();
            decimated.ShallowCopy(loop.GetOutput());

            Debug.WriteLine("After decimation" + Environment.NewLine + "------------");

            Debug.WriteLine("There are " + decimated.GetNumberOfPoints() + " points.");
            Debug.WriteLine("There are " + decimated.GetNumberOfPolys() + " polygons.");

            vtkPolyDataMapper inputMapper = vtkPolyDataMapper.New();
#if VTK_MAJOR_VERSION_5
         inputMapper.SetInputConnection(input.GetProducerPort());
#else
            inputMapper.SetInput(input);
#endif
            vtkActor inputActor = vtkActor.New();
            inputActor.SetMapper(inputMapper);

            vtkPolyDataMapper decimatedMapper = vtkPolyDataMapper.New();
#if VTK_MAJOR_VERSION_5
         decimatedMapper.SetInputConnection(decimated.GetProducerPort());
#else
            decimatedMapper.SetInput(decimated);
#endif
            vtkActor decimatedActor = vtkActor.New();
            decimatedActor.SetMapper(decimatedMapper);

            vtkRenderWindow renderWindow = myRenderWindowControl.RenderWindow;
            this.Size = new System.Drawing.Size(612, 352);

            // Define viewport ranges
            // (xmin, ymin, xmax, ymax)
            double[] leftViewport = new double[] { 0.0, 0.0, 0.5, 1.0 };
            double[] rightViewport = new double[] { 0.5, 0.0, 1.0, 1.0 };

            // Setup both renderers
            vtkRenderer leftRenderer = vtkRenderer.New();
            renderWindow.AddRenderer(leftRenderer);
            leftRenderer.SetViewport(leftViewport[0], leftViewport[1], leftViewport[2], leftViewport[3]);
            leftRenderer.SetBackground(.6, .5, .4);

            vtkRenderer rightRenderer = vtkRenderer.New();
            renderWindow.AddRenderer(rightRenderer);
            rightRenderer.SetViewport(rightViewport[0], rightViewport[1], rightViewport[2], rightViewport[3]);
            rightRenderer.SetBackground(.4, .5, .6);

            // Add the sphere to the left and the cube to the right
            leftRenderer.AddActor(inputActor);
            rightRenderer.AddActor(decimatedActor);
            leftRenderer.ResetCamera();
            rightRenderer.ResetCamera();
            renderWindow.Render();
        }

        private void Decimation()
        {
            vtkSphereSource sphereSource = vtkSphereSource.New();
            sphereSource.Update();

            vtkPolyData input = vtkPolyData.New();
            input.ShallowCopy(sphereSource.GetOutput());

            Debug.WriteLine("Before decimation" + Environment.NewLine + "------------");
            Debug.WriteLine("There are " + input.GetNumberOfPoints() + " points.");
            Debug.WriteLine("There are " + input.GetNumberOfPolys() + " polygons.");

            vtkDecimatePro decimate = vtkDecimatePro.New();
#if VTK_MAJOR_VERSION_5
         decimate.SetInputConnection(input.GetProducerPort());
#else
            decimate.SetInput(input);
#endif
            //decimate.SetTargetReduction(.99); //99% reduction (if there was 100 triangles, now there will be 1)
            decimate.SetTargetReduction(.50); //10% reduction (if there was 100 triangles, now there will be 50)
            decimate.Update();

            vtkPolyData decimated = vtkPolyData.New();
            decimated.ShallowCopy(decimate.GetOutput());

            Debug.WriteLine("After decimation" + Environment.NewLine + "------------");

            Debug.WriteLine("There are " + decimated.GetNumberOfPoints() + " points.");
            Debug.WriteLine("There are " + decimated.GetNumberOfPolys() + " polygons.");

            vtkPolyDataMapper inputMapper = vtkPolyDataMapper.New();
#if VTK_MAJOR_VERSION_5
         inputMapper.SetInputConnection(input.GetProducerPort());
#else
            inputMapper.SetInput(input);
#endif
            vtkActor inputActor = vtkActor.New();
            inputActor.SetMapper(inputMapper);

            vtkPolyDataMapper decimatedMapper = vtkPolyDataMapper.New();
#if VTK_MAJOR_VERSION_5
         decimatedMapper.SetInputConnection(decimated.GetProducerPort());
#else
            decimatedMapper.SetInput(decimated);
#endif
            vtkActor decimatedActor = vtkActor.New();
            decimatedActor.SetMapper(decimatedMapper);

            vtkRenderWindow renderWindow = myRenderWindowControl.RenderWindow;
            this.Size = new System.Drawing.Size(612, 352);

            // Define viewport ranges
            // (xmin, ymin, xmax, ymax)
            double[] leftViewport = new double[] { 0.0, 0.0, 0.5, 1.0 };
            double[] rightViewport = new double[] { 0.5, 0.0, 1.0, 1.0 };

            // Setup both renderers
            vtkRenderer leftRenderer = vtkRenderer.New();
            renderWindow.AddRenderer(leftRenderer);
            leftRenderer.SetViewport(leftViewport[0], leftViewport[1], leftViewport[2], leftViewport[3]);
            leftRenderer.SetBackground(.6, .5, .4);

            vtkRenderer rightRenderer = vtkRenderer.New();
            renderWindow.AddRenderer(rightRenderer);
            rightRenderer.SetViewport(rightViewport[0], rightViewport[1], rightViewport[2], rightViewport[3]);
            rightRenderer.SetBackground(.4, .5, .6);

            // Add the sphere to the left and the cube to the right
            leftRenderer.AddActor(inputActor);
            rightRenderer.AddActor(decimatedActor);
            leftRenderer.ResetCamera();
            rightRenderer.ResetCamera();
            renderWindow.Render();
        }

        private void WeightedTransformFilter()
        {
            // Use a sphere as a basis of the shape
            vtkSphereSource sphere = vtkSphereSource.New();
            sphere.SetPhiResolution(40);
            sphere.SetThetaResolution(40);
            sphere.Update();

            vtkPolyData sphereData = sphere.GetOutput();

            // Create a data array to hold the weighting coefficients
            vtkFloatArray tfarray = vtkFloatArray.New();
            int npoints = (int)sphereData.GetNumberOfPoints();
            tfarray.SetNumberOfComponents(2);
            tfarray.SetNumberOfTuples(npoints);

            // Parameterize the sphere along the z axis, and fill the weights
            // with (1.0-a, a) to linearly interpolate across the shape
            IntPtr pPoint = Marshal.AllocHGlobal(Marshal.SizeOf(typeof(double)) * 3);
            double[] point = new double[3];
            for (int i = 0; i<npoints; i++)
            {
                sphereData.GetPoint(i, pPoint);
                Marshal.Copy(pPoint, point, 0, 3);
                double x = point[0];
                double y = point[1];
                double z = point[2];

                double zn = z + 0.5;
                double zn1 = 1.0 - zn;
                if (zn<1.0)
                    zn = 1.0;
                if (zn1<0.0)
                    zn1 = 0.0;

                tfarray.SetComponent(i, 0, zn1);
                tfarray.SetComponent(i, 1, zn);
            }
            Marshal.FreeHGlobal(pPoint);

            // Create field data to hold the array, and bind it to the sphere
            vtkFieldData fd = vtkFieldData.New();
            tfarray.SetName("weights");
            sphereData.GetPointData().AddArray(tfarray);

            // Use an ordinary transform to stretch the shape
            vtkTransform stretch = vtkTransform.New();
            stretch.Scale(1, 1, 3.2);

            vtkTransformFilter stretchFilter = vtkTransformFilter.New();
            stretchFilter.SetInputConnection(sphereData.GetProducerPort());
            stretchFilter.SetTransform(stretch);

            // Now, for the weighted transform stuff
            vtkWeightedTransformFilter weightedTrans = vtkWeightedTransformFilter.New();

            // Create two transforms to interpolate between
            vtkTransform identity = vtkTransform.New();
            identity.Identity();

            vtkTransform rotated = vtkTransform.New();
            double rotatedAngle = 45;
            rotated.RotateX(rotatedAngle);

            weightedTrans.SetNumberOfTransforms(2);
            weightedTrans.SetTransform(identity, 0);
            weightedTrans.SetTransform(rotated, 1);
            // which data array should the filter use ?
            weightedTrans.SetWeightArray("weights");

            weightedTrans.SetInputConnection(stretchFilter.GetOutputPort());

            vtkPolyDataMapper weightedTransMapper = vtkPolyDataMapper.New();
            weightedTransMapper.SetInputConnection(weightedTrans.GetOutputPort());
            vtkActor weightedTransActor = vtkActor.New();
            weightedTransActor.SetMapper(weightedTransMapper);
            weightedTransActor.GetProperty().SetDiffuseColor(0.8, 0.8, 0.1);
            weightedTransActor.GetProperty().SetRepresentationToSurface();

            // get a reference to the renderwindow of our renderWindowControl1
            vtkRenderWindow renderWindow = myRenderWindowControl.RenderWindow;
            // renderer
            vtkRenderer renderer = renderWindow.GetRenderers().GetFirstRenderer();
            // set background color
            renderer.SetBackground(0.2, 0.3, 0.4);
            // add our actor to the renderer
            renderer.AddActor(weightedTransActor);

            renderer.ResetCamera();
            renderer.GetActiveCamera().Azimuth(90);
            renderer.GetActiveCamera().Dolly(1);
        }

        private void vtkPolyDataConnectivityFilter_LargestRegion()
        {
            // Small sphere
            vtkSphereSource sphereSource1 = vtkSphereSource.New();
            sphereSource1.Update();

            // Large sphere
            vtkSphereSource sphereSource2 = vtkSphereSource.New();
            sphereSource2.SetRadius(10);
            sphereSource2.SetCenter(25, 0, 0);
            sphereSource2.SetThetaResolution(10);
            sphereSource2.SetPhiResolution(10);
            sphereSource2.Update();

            vtkAppendPolyData appendFilter = vtkAppendPolyData.New();
            appendFilter.AddInputConnection(sphereSource1.GetOutputPort());
            appendFilter.AddInputConnection(sphereSource2.GetOutputPort());
            appendFilter.Update();

            vtkPolyDataConnectivityFilter connectivityFilter = vtkPolyDataConnectivityFilter.New();
            connectivityFilter.SetInputConnection(appendFilter.GetOutputPort());
            connectivityFilter.SetExtractionModeToLargestRegion();
            connectivityFilter.Update();

            // Create a mapper and actor for original data
            vtkPolyDataMapper originalMapper = vtkPolyDataMapper.New();
            originalMapper.SetInputConnection(appendFilter.GetOutputPort());
            originalMapper.Update();

            vtkActor originalActor = vtkActor.New();
            originalActor.SetMapper(originalMapper);

            // Create a mapper and actor for extracted data
            vtkPolyDataMapper extractedMapper = vtkPolyDataMapper.New();
            extractedMapper.SetInputConnection(connectivityFilter.GetOutputPort());
            extractedMapper.Update();

            vtkActor extractedActor = vtkActor.New();
            extractedActor.GetProperty().SetColor(1, 0, 0);
            extractedActor.SetMapper(extractedMapper);
            // get a reference to the renderwindow of our renderWindowControl1
            vtkRenderWindow renderWindow = myRenderWindowControl.RenderWindow;
            // renderer
            vtkRenderer renderer = renderWindow.GetRenderers().GetFirstRenderer();
            // set background color
            renderer.SetBackground(0.2, 0.3, 0.4);
            // add our actor to the renderer
            renderer.AddActor(originalActor);
            renderer.AddActor(extractedActor);
        }


        private void PolygonalSurfaceContourLineInterpolator()
        {
            vtkPolyData polyData;
            vtkSphereSource sphereSource = vtkSphereSource.New();
            sphereSource.SetThetaResolution(40);
            sphereSource.SetPhiResolution(20);
            sphereSource.Update();

            polyData = sphereSource.GetOutput();
            // The Dijkstra interpolator will not accept cells that aren't triangles
            vtkTriangleFilter triangleFilter = vtkTriangleFilter.New();
#if VTK_MAJOR_VERSION_5
         triangleFilter.SetInput( polyData );
#else
            triangleFilter.SetInput(polyData);
#endif
            triangleFilter.Update();

            vtkPolyData pd = triangleFilter.GetOutput();

            //Create a mapper and actor
            vtkPolyDataMapper mapper = vtkPolyDataMapper.New();
            mapper.SetInputConnection(triangleFilter.GetOutputPort());

            vtkActor actor = vtkActor.New();
            actor.SetMapper(mapper);
            actor.GetProperty().SetInterpolationToFlat();

            // get a reference to the renderwindow of our renderWindowControl1
            vtkRenderWindow renderWindow = myRenderWindowControl.RenderWindow;
            // renderer
            vtkRenderer renderer = renderWindow.GetRenderers().GetFirstRenderer();
            // set background color
            renderer.SetBackground(0.3, 0.4, 0.5);
            // add our actor to the renderer
            renderer.AddActor(actor);

            // Here comes the contour widget stuff.....
            vtkContourWidget contourWidget = vtkContourWidget.New();
            contourWidget.SetInteractor(renderWindow.GetInteractor());
            vtkOrientedGlyphContourRepresentation rep =
               vtkOrientedGlyphContourRepresentation.SafeDownCast(
                  contourWidget.GetRepresentation());
            rep.GetLinesProperty().SetColor(1, 0.2, 0);
            rep.GetLinesProperty().SetLineWidth(3.0f);

            vtkPolygonalSurfacePointPlacer pointPlacer =
               vtkPolygonalSurfacePointPlacer.New();
            pointPlacer.AddProp(actor);
            pointPlacer.GetPolys().AddItem(pd);
            rep.SetPointPlacer(pointPlacer);

            vtkPolygonalSurfaceContourLineInterpolator interpolator =
               vtkPolygonalSurfaceContourLineInterpolator.New();
            interpolator.GetPolys().AddItem(pd);
            rep.SetLineInterpolator(interpolator);

            renderWindow.Render();
            contourWidget.EnabledOn();
        }

        private void DrawPoint()
        {
            // Create the geometry of the points (the coordinate)
            vtkPoints points = vtkPoints.New();
            double[,] p = new double[,] 
            {
                {1.0, 2.0, 3.0},
                {3.0, 1.0, 2.0},
                {2.0, 3.0, 1.0},
                {1.0, 3.0, 3.0}
            };

            // Create topology of the points (a vertex per point)
            vtkCellArray vertices = vtkCellArray.New();
            int nPts = 4;

            int[] ids = new int[nPts];
             for(int i = 0; i<nPts; i++)
                ids[i] = (int)points.InsertNextPoint(p[i, 0], p[i, 1], p[i, 2]);

             int size = Marshal.SizeOf(typeof(int)) * nPts;
            IntPtr pIds = Marshal.AllocHGlobal(size);
            Marshal.Copy(ids, 0, pIds, nPts);
             vertices.InsertNextCell(nPts, pIds);
             Marshal.FreeHGlobal(pIds);

             // Create a polydata object
             vtkPolyData pointPoly = vtkPolyData.New();

            // Set the points and vertices we created as the geometry and topology of the polydata
            pointPoly.SetPoints(points);
             pointPoly.SetVerts(vertices);

             // Visualize
             vtkPolyDataMapper mapper = vtkPolyDataMapper.New();
             mapper.SetInput(pointPoly);

            vtkActor actor = vtkActor.New();
            actor.SetMapper(mapper);
            actor.GetProperty().SetPointSize(10);
            vtkRenderWindow renderWindow = myRenderWindowControl.RenderWindow;
            vtkRenderer renderer = renderWindow.GetRenderers().GetFirstRenderer();
            renderer.SetBackground(0.3, 0.2, 0.1);
            renderer.AddActor(actor);
        }

        private void DrawTriangle()
        {
            //创建点数据
            vtkPoints points = vtkPoints.New();
            points.InsertNextPoint(1.0, 0.0, 0.0);
            points.InsertNextPoint(0.0, 1.0, 0.0);
            points.InsertNextPoint(0.0, 0.0, 0.0);

            //每两个坐标之间分别创建一条线
            //SetId()的第一个参数是线段的端点ID，第二参数是连接的的点的ID
            vtkLine line0 = vtkLine.New();
            line0.GetPointIds().SetId(0, 0);
            line0.GetPointIds().SetId(1, 1);

            vtkLine line1 = vtkLine.New();
            line1.GetPointIds().SetId(0, 1);
            line1.GetPointIds().SetId(1, 2);

            vtkLine line2 = vtkLine.New();
            line2.GetPointIds().SetId(0, 2);
            line2.GetPointIds().SetId(1, 0);

            //创建单元数组，用于存储以上创建的线段
            vtkCellArray lines = vtkCellArray.New();
            lines.InsertNextCell(line0);
            lines.InsertNextCell(line1);
            lines.InsertNextCell(line2);

            //将点和线加入数据集中，前者定义数据集的几何结构，后者定义拓扑结构
            //创建vtkPolyData类型的数据，是一种数据集
            vtkPolyData polyData = vtkPolyData.New();

            //将创建的点数据加入vtkPolyData数据里
            polyData.SetPoints(points);  //点数据定义了polydata数据集的几何结构。
            polyData.SetLines(lines);   //定义拓扑结构

            //显示数据
            vtkPolyDataMapper mapper = vtkPolyDataMapper.New();
            mapper.SetInput(polyData);
            vtkActor actor = vtkActor.New();
            actor.SetMapper(mapper);
            actor.GetProperty().SetColor(1.0, 0.0, 0.0);           
            vtkRenderWindow renWin = myRenderWindowControl.RenderWindow;
            vtkRenderer renderer = renWin.GetRenderers().GetFirstRenderer();
            renderer.SetBackground(1.0, 1.0, 1.0);
            renderer.AddActor(actor);


        }

        private void DrawCylinder()
        {
            vtkCylinderSource cylinderSource = vtkCylinderSource.New();
            //cylinderSource.SetHeight(3.0);
            //cylinderSource.SetRadius(0.1);
            cylinderSource.SetResolution(8);
            vtkPolyDataMapper cylinderMapper = vtkPolyDataMapper.New();
            cylinderMapper.SetInputConnection(cylinderSource.GetOutputPort());
            vtkActor cylinderActor = vtkActor.New();
            cylinderActor.SetMapper(cylinderMapper);
            //cylinderActor.GetProperty().SetColor(1.0, 0.0, 0.0);
            cylinderActor.GetProperty().SetColor((float)Color.Tomato.R/256, (float)Color.Tomato.G/256, (float)Color.Tomato.B/256);
            cylinderActor.RotateX(30.0);
            cylinderActor.RotateY(-45.0);
            vtkRenderWindow renWin = myRenderWindowControl.RenderWindow;
            vtkRenderer renderer = renWin.GetRenderers().GetFirstRenderer();
            renderer.SetBackground(1.0, 1.0, 1.0);
            renderer.AddActor(cylinderActor);
        }

        private void DrawSphere(double radius)
        {
            vtkSphereSource sphereSource = vtkSphereSource.New();
            sphereSource.SetThetaResolution(8);
            sphereSource.SetPhiResolution(16);
            sphereSource.SetRadius(radius);

            vtkShrinkPolyData shrink = vtkShrinkPolyData.New();
            shrink.SetInputConnection(sphereSource.GetOutputPort());
            shrink.SetShrinkFactor(0.9);

            vtkPolyDataMapper sphereMapper = vtkPolyDataMapper.New();
            //sphereMapper.SetInputConnection(sphereSource.GetOutputPort());
            sphereMapper.SetInputConnection(shrink.GetOutputPort());

            vtkActor sphereActor = vtkActor.New();
            sphereActor.SetMapper(sphereMapper);
            sphereActor.GetProperty().SetColor(1, 0, 0);

            vtkRenderer sphereRender = vtkRenderer.New();
            vtkRenderWindow renWin = myRenderWindowControl.RenderWindow;
            renWin.AddRenderer(sphereRender);

            sphereRender.AddActor(sphereActor);
            sphereRender.SetBackground(0.0, 0.0, 1.0);
        }

        private void ReadStlFileDraw()
        {
            vtkSTLReader stlReader = vtkSTLReader.New();
            //stlReader.SetFileName(@"..\..\Data\anchor_hook.stl");
            stlReader.SetFileName(@"..\..\Data\42400-IDGH.stl");
            vtkShrinkPolyData shrink = vtkShrinkPolyData.New();
            shrink.SetInputConnection(stlReader.GetOutputPort());
            shrink.SetShrinkFactor(0.85);

            vtkPolyDataMapper stlMapper = vtkPolyDataMapper.New();
            //stlMapper.SetInputConnection(stlReader.GetOutputPort());
            stlMapper.SetInputConnection(shrink.GetOutputPort());

            vtkLODActor stlActor = vtkLODActor.New();
            stlActor.SetMapper(stlMapper);
            stlActor.GetProperty().SetColor((float)Color.Orange.R/256, (float)Color.Orange.G / 256,(float)Color.Orange.B / 256);  //设置actor的颜色

            vtkRenderer stlRender = vtkRenderer.New();
            vtkRenderWindow renWin = myRenderWindowControl.RenderWindow;
            renWin.AddRenderer(stlRender);

            stlRender.AddActor(stlActor);
            stlRender.SetBackground(1.0, 1.0, 1.0);

        }

        private void DrawAssembly()
        {
            //Create four parts: a top level assembly and three primitives

            vtkSphereSource sphereSource = vtkSphereSource.New();
            vtkPolyDataMapper sphereMapper = vtkPolyDataMapper.New();
            sphereMapper.SetInputConnection(sphereSource.GetOutputPort());
            vtkActor sphereActor = vtkActor.New();
            sphereActor.SetMapper(sphereMapper);
            sphereActor.SetOrigin(2, 1, 3);
            sphereActor.RotateY(6);
            sphereActor.SetPosition(2.25, 0, 0);
            sphereActor.GetProperty().SetColor(1, 0, 1);

            vtkCubeSource cubeSource = vtkCubeSource.New();
            vtkPolyDataMapper cubeMapper = vtkPolyDataMapper.New();
            cubeMapper.SetInputConnection(cubeSource.GetOutputPort());
            vtkActor cubeActor = vtkActor.New();
            cubeActor.SetMapper(cubeMapper);
            cubeActor.SetPosition(0, 2.25, 0);
            cubeActor.GetProperty().SetColor(0, 0, 1);

            vtkConeSource coneSource = vtkConeSource.New();
            vtkPolyDataMapper coneMapper = vtkPolyDataMapper.New();
            coneMapper.SetInputConnection(coneSource.GetOutputPort());
            vtkActor coneActor = vtkActor.New();
            coneActor.SetMapper(coneMapper);
            coneActor.SetPosition(0, 0, 2.25);
            coneActor.GetProperty().SetColor(0, 1, 0);

            vtkCylinderSource cylinderSource = vtkCylinderSource.New();
            vtkPolyDataMapper cylinderMapper = vtkPolyDataMapper.New();
            cylinderMapper.SetInputConnection(cylinderSource.GetOutputPort());
            vtkActor cylinderActor = vtkActor.New();
            cylinderActor.SetMapper(cylinderMapper);
            //cylinderActor.SetPosition(0, 0, 0);
            cylinderActor.GetProperty().SetColor(1, 0, 0);

            vtkAssembly assembly = vtkAssembly.New();
            assembly.AddPart(cylinderActor);
            assembly.AddPart(sphereActor);
            assembly.AddPart(cubeActor);
            assembly.AddPart(coneActor);
            assembly.SetOrigin(5, 10, 5);
            assembly.AddPosition(5, 0, 0);
            assembly.RotateX(15);

            vtkRenderer renderer = vtkRenderer.New();
            vtkRenderWindow renWin = myRenderWindowControl.RenderWindow;

            renWin.AddRenderer(renderer);
            renderer.AddActor(assembly);
            renderer.AddActor(coneActor);

        }

        private void DrawTexturePlane()
        {
            //load in the texture map
            vtkBMPReader bmpReader = vtkBMPReader.New();           
            bmpReader.SetFileName(@"..\..\Data\masonry.bmp");
            vtkTexture atext = vtkTexture.New();
            atext.SetInputConnection(bmpReader.GetOutputPort());
            atext.InterpolateOn();

            //create a plane source and actor
            vtkPlaneSource plane = vtkPlaneSource.New();
            plane.SetPoint1(0, 0, 0);
            vtkPolyDataMapper planeMapper = vtkPolyDataMapper.New();
            planeMapper.SetInputConnection(plane.GetOutputPort());
            vtkActor planeActor = vtkActor.New();
            planeActor.SetMapper(planeMapper);
            planeActor.SetTexture(atext);

            vtkRenderer renderer = vtkRenderer.New();
            vtkRenderWindow renWin = myRenderWindowControl.RenderWindow;

            renWin.AddRenderer(renderer);
            renderer.AddActor(planeActor);
        }

        private void DrawRainBow()
        {
            //# First create pipeline a simple pipeline that reads a structure grid
            //# and then extracts a plane from the grid. The plane will be colored
            //# differently by using different lookup tables.
            //#
            //# Note: the Update method is manually invoked because it causes the
            //# reader to read; later on we use the output of the reader to set
            //# a range for the scalar values.
            vtkMultiBlockPLOT3DReader pl3d = vtkMultiBlockPLOT3DReader.New();
            pl3d.SetXYZFileName(@"..\..\Data\combxyz.bin");
            pl3d.SetQFileName(@"..\..\Data\combq.bin");
            pl3d.SetScalarFunctionNumber(100);
            pl3d.SetVectorFunctionNumber(202);
            pl3d.Update();
            vtkDataObject pl3d_output = pl3d.GetOutput().GetBlock(0);

            vtkStructuredGridGeometryFilter planeFilter = vtkStructuredGridGeometryFilter.New();
            planeFilter.SetInput(pl3d_output);
            planeFilter.SetExtent(1, 100, 1, 100, 7, 7);
            vtkLookupTable lut = vtkLookupTable.New();
            vtkPolyDataMapper planeMapper = vtkPolyDataMapper.New();
            planeMapper.SetLookupTable(lut);
            planeMapper.SetInputConnection(planeFilter.GetOutputPort());
            //planeMapper.SetScalarRange(pl3d_output.)
            vtkActor planeActor = vtkActor.New();
            planeActor.SetMapper(planeMapper);

            //this creates an outline around the data
            vtkStructuredGridOutlineFilter outlineFilter = vtkStructuredGridOutlineFilter.New();
            outlineFilter.SetInput(pl3d_output);
            vtkPolyDataMapper outlineMapper = vtkPolyDataMapper.New();
            outlineMapper.SetInputConnection(outlineFilter.GetOutputPort());
            vtkActor outlineActor = vtkActor.New();
            outlineActor.SetMapper(outlineMapper);

            //Much of the following is commented out. To try different lookup tables.
            //This create a black to white lut
            //lut.SetHueRange(0, 0);
            //lut.SetSaturationRange(0, 0);
            //lut.SetValueRange(0.2, 1.0);

            //This creates a red to blue lut
            //lut.SetHueRange(0.0, 0.677);

            //This creates a blue to red lue
            lut.SetHueRange(0.667, 0.0);

            //This creates a weird effect. the Build() method cause lookup
            //table to allocate memory and create a table based on the correct
            //hue, saturatioin, value, and alpha range. Here we then 
            //manully overwrite the value generated by the Build() method.
            lut.SetNumberOfColors(256);
            lut.Build();
            for(int i=0;i<16;i++)
            {
                lut.SetTableValue(i * 16, (float)Color.Red.R / 256, (float)Color.Red.G / 256, (float)Color.Red.B / 256, 1);
                lut.SetTableValue(i * 16+1, (float)Color.Green.R / 256, (float)Color.Green.G / 256, (float)Color.Green.B / 256, 1);
                lut.SetTableValue(i * 16+2, (float)Color.Blue.R / 256, (float)Color.Blue.G / 256, (float)Color.Blue.B / 256, 1);
                lut.SetTableValue(i * 16+3, (float)Color.Black.R / 256, (float)Color.Black.G / 256, (float)Color.Black.B / 256, 1);
            }


            //Create the renderwindow, the render and both actors
            vtkRenderer ren = vtkRenderer.New();
            vtkRenderWindow renWin = myRenderWindowControl.RenderWindow;
            renWin.AddRenderer(ren);

            //Add the actors to the renderer, set the backgroud
            ren.AddActor(outlineActor);
            ren.AddActor(planeActor);

            ren.SetBackground(0.1, 0.2, 0.4);
            ren.TwoSidedLightingOff();
        }
        
        private void DrawVisQuad()
        {
            //# This example demonstrates the use of the contour filter, and the use of
            //# the vtkSampleFunction to generate a volume of data samples from an
            //# implicit function.

            //# VTK supports implicit functions of the form f(x,y,z)=constant. These
            //# functions can represent things spheres, cones, etc. Here we use a
            //# general form for a quadric to create an elliptical data field.

            vtkQuadric quadricFunction = vtkQuadric.New();
            quadricFunction.SetCoefficients(0.5, 1, 0.2, 0, 0.1, 0, 0, 0.2, 0, 0);

            //# vtkSampleFunction samples an implicit function over the x-y-z range
            //# specified (here it defaults to -1,1 in the x,y,z directions).
            vtkSampleFunction sample = vtkSampleFunction.New();
            sample.SetSampleDimensions(30, 30, 30);
            sample.SetImplicitFunction(quadricFunction);

            //# Create five surfaces F(x,y,z) = constant between range specified. The
            //# GenerateValues() method creates n isocontour values between the range
            //# specified.
            vtkContourFilter contourFilter = vtkContourFilter.New();
            contourFilter.SetInputConnection(sample.GetOutputPort());
            contourFilter.GenerateValues(10, 0.0, 1.2);

            vtkPolyDataMapper contMapper = vtkPolyDataMapper.New();
            contMapper.SetInputConnection(contourFilter.GetOutputPort());
            contMapper.SetScalarRange(0.0, 1.2);

            vtkActor conActor = vtkActor.New();
            conActor.SetMapper(contMapper);

            //We'll put a simple outline around the data
            vtkOutlineFilter outline = vtkOutlineFilter.New();
            outline.SetInputConnection(sample.GetOutputPort());

            vtkPolyDataMapper outlineMapper = vtkPolyDataMapper.New();
            outlineMapper.SetInputConnection(outline.GetOutputPort());

            vtkActor outlineActor = vtkActor.New();
            outlineActor.SetMapper(outlineMapper);
            outlineActor.GetProperty().SetColor(0, 0, 0);

            //The usual rendering stuff
            vtkRenderer ren = vtkRenderer.New();
            vtkRenderWindow renWin = myRenderWindowControl.RenderWindow;
            //vtkRenderWindow renWin = vtkRenderWindow.New();
            renWin.AddRenderer(ren);
            //vtkRenderWindowInteractor iren = vtkRenderWindowInteractor.New();
            //iren.SetRenderWindow(renWin);

            ren.SetBackground(1, 1, 1);
            ren.AddActor(conActor);
            ren.AddActor(outlineActor);

            //iren.Initialize();
            //renWin.Render();
            //iren.Start();

        }

        private void DrawBuildUGrid()
        {
            //# This example shows how to manually construct unstructured grids
            //# using C#.  Unstructured grids require explicit point and cell
            //# representations, so every point and cell must be created, and then
            //# added to the vtkUnstructuredGrid instance.

            //# Create several unstructured grids each containing a cell of a
            //# different type.

            //create voxel
            vtkPoints voxelPoints = vtkPoints.New();
            voxelPoints.SetNumberOfPoints(8);
            voxelPoints.InsertPoint(0, 0, 0, 0);
            voxelPoints.InsertPoint(1, 1, 0, 0);
            voxelPoints.InsertPoint(2, 0, 1, 0);
            voxelPoints.InsertPoint(3, 1, 1, 0);
            voxelPoints.InsertPoint(4, 0, 0, 1);
            voxelPoints.InsertPoint(5, 1, 0, 1);
            voxelPoints.InsertPoint(6, 0, 1, 1);
            voxelPoints.InsertPoint(7, 1, 1, 1);
            vtkVoxel aVoxel = vtkVoxel.New();
            aVoxel.GetPointIds().SetId(0, 0);
            aVoxel.GetPointIds().SetId(1, 1);
            aVoxel.GetPointIds().SetId(2, 2);
            aVoxel.GetPointIds().SetId(3, 3);
            aVoxel.GetPointIds().SetId(4, 4);
            aVoxel.GetPointIds().SetId(5, 5);
            aVoxel.GetPointIds().SetId(6, 6);
            aVoxel.GetPointIds().SetId(7, 7);
            vtkUnstructuredGrid aVoxelGrid = vtkUnstructuredGrid.New();
            aVoxelGrid.Allocate(1, 1);
            aVoxelGrid.InsertNextCell(aVoxel.GetCellType(), aVoxel.GetPointIds());
            aVoxelGrid.SetPoints(voxelPoints);
            vtkDataSetMapper aVoxelMapper = vtkDataSetMapper.New();
            aVoxelMapper.SetInput(aVoxelGrid);
            vtkActor aVoxelActor = vtkActor.New();
            aVoxelActor.SetMapper(aVoxelMapper);
            aVoxelActor.GetProperty().SetDiffuseColor(1, 0, 0);

            //create Hexahedron
            vtkPoints hexahedronPoints = new vtkPoints();
            hexahedronPoints.SetNumberOfPoints(8);
            hexahedronPoints.InsertPoint(0, 0, 0, 0);
            hexahedronPoints.InsertPoint(1, 1, 0, 0);
            hexahedronPoints.InsertPoint(2, 1, 1, 0);
            hexahedronPoints.InsertPoint(3, 0, 1, 0);
            hexahedronPoints.InsertPoint(4, 0, 0, 1);
            hexahedronPoints.InsertPoint(5, 1, 0, 1);
            hexahedronPoints.InsertPoint(6, 1, 1, 1);
            hexahedronPoints.InsertPoint(7, 0, 1, 1);
            vtkHexahedron aHexahedron = new vtkHexahedron();
            aHexahedron.GetPointIds().SetId(0, 0);
            aHexahedron.GetPointIds().SetId(1, 1);
            aHexahedron.GetPointIds().SetId(2, 2);
            aHexahedron.GetPointIds().SetId(3, 3);
            aHexahedron.GetPointIds().SetId(4, 4);
            aHexahedron.GetPointIds().SetId(5, 5);
            aHexahedron.GetPointIds().SetId(6, 6);
            aHexahedron.GetPointIds().SetId(7, 7);
            vtkUnstructuredGrid aHexahedronGrid = new vtkUnstructuredGrid();
            aHexahedronGrid.Allocate(1, 1);
            aHexahedronGrid.InsertNextCell(aHexahedron.GetCellType(), aHexahedron.GetPointIds());
            aHexahedronGrid.SetPoints(hexahedronPoints);

            vtkDataSetMapper aHexahedronMapper = new vtkDataSetMapper();
            aHexahedronMapper.SetInput(aHexahedronGrid);
            vtkActor aHexahedronActor = new vtkActor();
            aHexahedronActor.SetMapper(aHexahedronMapper);
            aHexahedronActor.AddPosition(2, 0, 0);
            aHexahedronActor.GetProperty().SetDiffuseColor(1, 1, 0);


            vtkRenderer render = vtkRenderer.New();
            vtkRenderWindow renWin = myRenderWindowControl.RenderWindow;
            renWin.AddRenderer(render);

            render.SetBackground(0, 0, 1);
            render.AddActor(aVoxelActor);
            render.AddActor(aHexahedronActor);

        }

        private void DrawTest()
        {
            vtkProp3D prop3D;
            vtkActor actor = vtkActor.New();
            vtkActor2D actor2D = vtkActor2D.New();
            vtkLODActor lODActor = vtkLODActor.New();
            vtkLODProp3D lodProp3d = vtkLODProp3D.New();
            vtkCamera camera = vtkCamera.New();
            vtkCameraActor cameraActor = vtkCameraActor.New();
            vtkLight light = vtkLight.New();
            vtkLightActor lightActor = vtkLightActor.New();
            vtkPicker picker = vtkPicker.New();
            vtkPointPicker pointPicker = vtkPointPicker.New();
            vtkCellPicker cellPicker = vtkCellPicker.New();
            vtkAreaPicker areaPicker = vtkAreaPicker.New();

            vtkAssembly assembly = vtkAssembly.New();
            vtkConeSource coneSource = vtkConeSource.New();
            vtkCone cone = vtkCone.New();

            vtkArcSource arcSource = vtkArcSource.New();
            vtkLineSource lineSource = vtkLineSource.New();
            vtkPointSource pointSource = vtkPointSource.New();

            vtkPolyData polyData = vtkPolyData.New();
            vtkArrayReader arrayReader = vtkArrayReader.New();
            //vtkArrayDataReader arrayDataReader = vtkArrayDataReader.New();
            vtkArrayWriter arrayWriter = vtkArrayWriter.New();
            vtkRenderWindowInteractor renderWindowInteractor = vtkRenderWindowInteractor.New();
            //vtkRenderWindowInteractor3D renderWindowInteractor3D = vtkRenderWindowInteractor3D.New();
            vtkInteractorStyle interactorStyle = vtkInteractorStyle.New();
            //vtkInteractorStyle3D interactorStyle3D = vtkInteractorStyle3D.New();
            vtkInteractorStyleFlight interactorStyleFlight = vtkInteractorStyleFlight.New();
            vtkInteractorStyleTrackball interactorStyleTrackball = vtkInteractorStyleTrackball.New();

            vtkVolume volume = vtkVolume.New();
            vtkVolumeMapper volumeMapper;
            vtkSmartVolumeMapper smartVolumeMapper = vtkSmartVolumeMapper.New();
            vtkUnstructuredGridVolumeMapper unstructuredGridVolumeMapper;
            vtkUnstructuredGridVolumeRayCastMapper unstructuredGridVolumeRayCastMapper = vtkUnstructuredGridVolumeRayCastMapper.New();
            vtkGPUVolumeRayCastMapper gPUVolumeRayCastMapper = vtkGPUVolumeRayCastMapper.New();
            vtkVolumeRayCastMapper volumeRayCastMapper = vtkVolumeRayCastMapper.New();
            vtkFixedPointVolumeRayCastMapper pointVolumeRayCastMapper = vtkFixedPointVolumeRayCastMapper.New();
            vtkOpenGLGPUVolumeRayCastMapper openGLGPUVolumeRayCastMapper = vtkOpenGLGPUVolumeRayCastMapper.New();
            vtkVolumeProperty volumeProperty = vtkVolumeProperty.New();

            vtkTexture texture = vtkTexture.New();
            vtkCoordinate coordinate = vtkCoordinate.New();
            vtkImageData vtkImage = vtkImageData.New();

            vtkBMPReader bMPReader = vtkBMPReader.New();
            vtkJPEGReader jPEGReader = vtkJPEGReader.New();
            vtkPNGReader pNGReader = vtkPNGReader.New();
            vtkTIFFReader tIFFReader = vtkTIFFReader.New();
            vtkOBJReader oBJReader = vtkOBJReader.New();


            vtkContourFilter contourFilter = vtkContourFilter.New();
            vtkSynchronizedTemplates2D synchronizedTemplates2D = vtkSynchronizedTemplates2D.New();
            vtkSynchronizedTemplates3D synchronizedTemplates3D = vtkSynchronizedTemplates3D.New();
            vtkSynchronizedTemplatesCutter3D synchronizedTemplatesCutter3D = vtkSynchronizedTemplatesCutter3D.New();

            vtkImageMapper imageMapper = vtkImageMapper.New();
            //vtkImageSliceMapper imageSliceMapper = vtkImageSliceMapper.New();
            //vtkImageResliceMapper imageResliceMapper = vtkImageResliceMapper.New();

            vtkStructuredGridReader structuredGridReader = vtkStructuredGridReader.New(); 
            vtkRungeKutta4 integ = vtkRungeKutta4.New();
            vtkStreamTracer streamer = vtkStreamTracer.New();
            vtkTubeFilter streamTube = vtkTubeFilter.New();
            vtkRuledSurfaceFilter ruledSurfaceFilter = vtkRuledSurfaceFilter.New();
            vtkPlane plane = vtkPlane.New();
            vtkCutter cutter = new vtkCutter();
            vtkMergeFilter mergeFilter = vtkMergeFilter.New();
            vtkImageLuminance imageLuminance = vtkImageLuminance.New();
            vtkImageDataGeometryFilter imageDataGeometryFilter = vtkImageDataGeometryFilter.New();
            vtkWarpScalar warpScalar = vtkWarpScalar.New();
            vtkWarpVector warpVector = vtkWarpVector.New();

        }

        private void panel1_Resize(object sender, EventArgs e)
        {
            myRenderWindowControl.SetBounds(0, 0, panel1.Width, panel1.Height);
            
        }
    }
}
