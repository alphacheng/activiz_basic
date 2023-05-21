using Kitware.VTK;
using System;

namespace point_vs_cell
{
    internal class Program
    {
        /*
         
2 数据前提
** 必备数据：**离散点数据（X1，Y1，Z1）、（X2，Y2，Z2）、（X3，Y3，Z3）…
**中间数据：**离散点三角剖分拓扑关系（即哪几个点构成了哪一个三角形），常用数据组成为（三角形序号，点1序号，点2序号，点3序号）

离散点数据就是我们最最基础的源数据了，这个是不可或缺的，也就是作为已知数据。对于中间数据（三角剖分拓扑关系数据）如果不是已知数据，我们有两种方案进行生成：
（1）基于Triangle工具进行三角剖分生成 （我使用的Triangle工具）
（2）基于VTK自身的Delaunay多边形过滤器（vtkDelaunay2D/vtkDelaunay3D）对点集（VtkPoints）进行过滤生成
现阶段中间数据作为已知数据进行处理。

3 构建思路
针对于我们有离散点和离散点三角剖分拓扑关系的数据前提去构建等值线和标注，有两种构建模式：
（1）构建点集——构建Cell单元集（拓扑关系）——构建结点标量集——构建模型
（2）构建点集——构建Cell单元集（拓扑关系）——构建Cell标量集——构建模型
第一种方案是将每一个结点看作一个标量单位，适用于对点元模型的渲染，即一个区域由一堆点元控制，每个点元（结点）有一个数据属性，例如一口井的水位值、一座烟囱的排放量…
第二种则将每一个Cell单元看作一个标量单位，适用于面元模型的渲染，即一片区域被划分为了一个个的单元且每一个单元有一个数据属性，例如每个单元的降雨量、事故发生率…
————————————————
版权声明：本文为CSDN博主「billy_gisboy」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
原文链接：https://blog.csdn.net/weixin_44188072/article/details/114983889
           
         */
        static void Main(string[] args)
        {

            /*构建网格结点*/
            vtkPoints points = new vtkPoints();
            points.InsertPoint(1, 0, 0, 0);
            points.InsertPoint(2, 1, 0, 0);
            points.InsertPoint(3, 2, 0, 0);
            points.InsertPoint(4, 0, 1, 0);
            points.InsertPoint(5, 1, 1, 0);
            points.InsertPoint(6, 2, 1, 0);
            points.InsertPoint(7, 0, 2, 0);
            points.InsertPoint(8, 1, 2, 0);
            points.InsertPoint(9, 2, 2, 0);
            points.InsertPoint(10, 3, 2, 0);
            points.InsertPoint(11, 3, 1, 0);

            /*构建网格单元*/
            vtkIdList idList1 = new vtkIdList();
            idList1.InsertNextId(1);
            idList1.InsertNextId(2);
            idList1.InsertNextId(4);
            vtkIdList idList2 = new vtkIdList();
            idList2.InsertNextId(4);
            idList2.InsertNextId(5);
            idList2.InsertNextId(7);
            vtkIdList idList3 = new vtkIdList();
            idList3.InsertNextId(2);
            idList3.InsertNextId(5);
            idList3.InsertNextId(4);
            vtkIdList idList4 = new vtkIdList();
            idList4.InsertNextId(3);
            idList4.InsertNextId(6);
            idList4.InsertNextId(5);
            vtkIdList idList5 = new vtkIdList();
            idList5.InsertNextId(2);
            idList5.InsertNextId(3);
            idList5.InsertNextId(5);
            vtkIdList idList6 = new vtkIdList();
            idList6.InsertNextId(6);
            idList6.InsertNextId(9);
            idList6.InsertNextId(8);
            vtkIdList idList7 = new vtkIdList();
            idList7.InsertNextId(5);
            idList7.InsertNextId(6);
            idList7.InsertNextId(8);
            vtkIdList idList8 = new vtkIdList();
            idList8.InsertNextId(5);
            idList8.InsertNextId(8);
            idList8.InsertNextId(7);
            vtkIdList idList9 = new vtkIdList();
            idList9.InsertNextId(9);
            idList9.InsertNextId(6);
            idList9.InsertNextId(11);
            vtkIdList idList10 = new vtkIdList();
            idList10.InsertNextId(9);
            idList10.InsertNextId(11);
            idList10.InsertNextId(10);

            vtkCellArray cellArray = new vtkCellArray();
            cellArray.InsertNextCell(idList1);
            cellArray.InsertNextCell(idList2);
            cellArray.InsertNextCell(idList3);
            cellArray.InsertNextCell(idList4);
            cellArray.InsertNextCell(idList5);
            cellArray.InsertNextCell(idList6);
            cellArray.InsertNextCell(idList7);
            cellArray.InsertNextCell(idList8);
            cellArray.InsertNextCell(idList9);
            cellArray.InsertNextCell(idList10);


            /*结点标量构建思路*/
            /*设置结点标量*/
            vtkFloatArray pointScalar = new vtkFloatArray();
            //这个插入方式必须和单元集（vtkCellArray）的插入方式一致，单元集插入Cell单元采用InsertNextCell自动递增插入
            pointScalar.InsertNextValue(1); //1
            pointScalar.InsertNextValue(2); //2
            pointScalar.InsertNextValue(3); //3
            pointScalar.InsertNextValue(2); //4
            pointScalar.InsertNextValue(3); //5
            pointScalar.InsertNextValue(4); //6
            pointScalar.InsertNextValue(3); //7
            pointScalar.InsertNextValue(4); //8
            pointScalar.InsertNextValue(5); //9
            pointScalar.InsertNextValue(4); //10
            pointScalar.InsertNextValue(5); //11

            double rangeMin = pointScalar.GetRange()[0];
            double rangeMax = pointScalar.GetRange()[1];

            /*构建表面几何模型*/
            //PolyData
            vtkPolyData surfacePolydata = vtkPolyData.New();
            surfacePolydata.SetPoints(points);
            surfacePolydata.SetPolys(cellArray); //设置数据类型为多边形，并设置多边形单元数据
            surfacePolydata.GetPointData().SetScalars(pointScalar); //设置点集标量
            //Actor：表面1 构建无规则网格单元的表面模型，效果是网格单元中的一个个面片，有颜色过渡
            vtkDataSetMapper surfaceMapper = new vtkDataSetMapper();
            surfaceMapper.SetInput(surfacePolydata);
            surfaceMapper.ScalarVisibilityOn(); //开启标量颜色渲染
            surfaceMapper.SetScalarRange(rangeMin, rangeMax);
            vtkLookupTable lookupTable = vtkLookupTable.New();
            lookupTable.SetRange(rangeMin, rangeMax);
            surfaceMapper.SetLookupTable(lookupTable); //此处可以根据标量范围设置自己需要的颜色映射表，其它的制图器Mapper也一样
            vtkActor surfaceActor = new vtkActor();
            surfaceActor.SetMapper(surfaceMapper);

            /*构建等值线*/
            //Fileter：等值线/面过滤器
            vtkContourFilter contourFilter = new vtkContourFilter();
            contourFilter.SetInput(surfacePolydata);
            contourFilter.GenerateValues(8, rangeMin, rangeMax);
            contourFilter.Update();
            //数据处理类：用于将离散的三角面片拼接为连续的等值面
            vtkStripper stripper = new vtkStripper();
            stripper.SetInput(contourFilter.GetOutput());
            stripper.Update();
            //Mapper
            vtkDataSetMapper contourMapper = new vtkDataSetMapper();
            contourMapper.SetInput(stripper.GetOutput());
            //Actor
            vtkActor isolineActor = new vtkActor();
            isolineActor.SetMapper(contourMapper);

            vtkScalarBarActor scalarBarActor = vtkScalarBarActor.New();
            scalarBarActor.SetLookupTable(lookupTable);


            /*等值线值标注*/
            //构建等值线标注几何数据
            //long linesCount = stripper.GetOutput().GetNumberOfLines(); //等值线的数目
            vtkPoints isolinePoints = stripper.GetOutput().GetPoints(); //等值线的点集
            vtkCellArray isolineCells = stripper.GetOutput().GetLines(); //等值线的线单元数组
            vtkDataArray isolineScalars = stripper.GetOutput().GetPointData().GetScalars(); //等值线的点标量数组
            //等值线标注Label数据
            vtkPolyData labelPolyData = new vtkPolyData(); //标注点几何数据
            vtkPoints labelPoints = new vtkPoints(); //标注点集
            vtkDoubleArray labelScalars = new vtkDoubleArray(); //标注点标量数组
            labelScalars.SetNumberOfComponents(1);
            //labelScalars.SetName("Isovalues");
            //构建等值线标注Label数据，此处将每一个等值线进行遍历，获取线上的一点作为等值线的标注点
            vtkIdList isoCell = new vtkIdList();
            while (isolineCells.GetNextCell(isoCell) > 0)
            {
                long idCount = isoCell.GetNumberOfIds(); //获取这条等值线单元的点个数
                //long samplePointIndex = (long)vtkMath.Random(0, idCount); //随机取该等值线上一个点
                long labelPointId = isoCell.GetId(idCount / 2); //获取该点的id作为标记点。此方法是根据点数量顺序排序序号Index获取，获取到这个cell中的指定序号Index的点的vtkId（个人认为）

                double[] midPoint = isolinePoints.GetPoint(labelPointId); //获取该点的坐标

                labelPoints.InsertNextPoint(midPoint[0], midPoint[1], midPoint[2]); //插入到标注点集
                labelScalars.InsertNextTuple1(isolineScalars.GetTuple1(labelPointId)); //取该点的标量值插入到标注点标量数组
            }
            labelPolyData.SetPoints(labelPoints); //设置标注点集
            labelPolyData.GetPointData().SetScalars(labelScalars); //设置点数据的标量数组
            labelPolyData.Update();
            //Mapper
            vtkLabeledDataMapper labelMapper = new vtkLabeledDataMapper();
            //labelMapper.SetFieldDataName("Isovalues");
            labelMapper.SetInput(labelPolyData);
            labelMapper.SetLabelModeToLabelScalars();
            labelMapper.SetLabelFormat("%6.2f");
            labelMapper.GetLabelTextProperty().SetColor(0, 0, 1);
            //Actor
            vtkActor2D isoLabelActor = new vtkActor2D();
            isoLabelActor.SetMapper(labelMapper);


            /*设置单元标量*/
            vtkFloatArray cellScalar = new vtkFloatArray();
            //这个插入方式必须和单元集（vtkCellArray）的插入方式一致，单元集插入cell单元采用InsertNextCell自动递增插入
            cellScalar.InsertNextValue(1); //1
            cellScalar.InsertNextValue(2); //2
            cellScalar.InsertNextValue(2); //3
            cellScalar.InsertNextValue(3); //4
            cellScalar.InsertNextValue(2); //5
            cellScalar.InsertNextValue(4); //6
            cellScalar.InsertNextValue(3); //7
            cellScalar.InsertNextValue(3); //8
            cellScalar.InsertNextValue(3); //9
            cellScalar.InsertNextValue(4); //10

            double rangeMin2 = cellScalar.GetRange()[0];
            double rangeMax2 = cellScalar.GetRange()[1];


            /*构建表面模型：效果是网格单元中的一个个面片，每个面片有这个cell标量对应的色值*/
            //Data：非结构化网格对象
            vtkUnstructuredGrid grid = new vtkUnstructuredGrid();
            grid.SetPoints(points);
            grid.SetCells(new vtkTriangle().GetCellType(), cellArray); //因为后续用到了三角单元数据处理类，处理的数据都是三角形，所以直接使用确定的cell类型
            grid.GetCellData().SetScalars(cellScalar); //设置标量
            //Mapper
            vtkDataSetMapper surfaceMapper2 = new vtkDataSetMapper();
            surfaceMapper2.SetInput(grid);
            //surfaceMapper.ScalarVisibilityOn();
            surfaceMapper2.SetScalarRange(rangeMin2, rangeMax2);
            //Actor
            vtkActor surfaceActor2 = new vtkActor();
            surfaceActor2.SetMapper(surfaceMapper2);


            ///vtkCellDataToPointData
            ///官方说明：
            ///（1）是一个用于转换单元格数据的过滤器
            ///（2）它将每个单元指定的数据转换为点数据（即在单元点指定的数据）
            ///（3）变换方法是基于使用特定点对所有单元的数据值求平均值，也可以选择将输入单元格数据传递到输出
            ///个人理解/问题：
            ///（1）为什么经过vtkCellDataToPointData处理后产生的模型颜色会均匀过渡？该工具类将grid的单元数据进行了处理，将标量在各个点进行了均匀的计算赋值（插值）

            /*构建表面模型：效果是网格单元中的一个个面片，但通过cell转点，实现了颜色的均匀过渡*/
            //Data:将构建的非结构化网格对象的每个单元转换为结点
            vtkCellDataToPointData cellDataToPointData = new vtkCellDataToPointData();
            cellDataToPointData.SetInput(grid);
            cellDataToPointData.PassCellDataOn();
            cellDataToPointData.Update();
            System.Console.WriteLine(  cellDataToPointData.GetOutput().GetScalarRange().ToString());

            //Mapper
            vtkDataSetMapper surfaceMapper3 = new vtkDataSetMapper();
            surfaceMapper3.SetInput(cellDataToPointData.GetOutput());
            //surfaceMapper2.SetScalarModeToUsePointData();
            //surfaceMapper2.ScalarVisibilityOn();
            surfaceMapper3.SetScalarRange(rangeMin2, rangeMax2);
            //Actor
            vtkActor surfaceActor3 = new vtkActor();
            surfaceActor3.SetMapper(surfaceMapper3);


            /*构建等值线*/
            //等值线/面过滤器
            vtkContourFilter contourFilter2 = new vtkContourFilter();
            contourFilter2.SetInput(cellDataToPointData.GetOutput());
            contourFilter2.GenerateValues(8, rangeMin2, rangeMax2);
            contourFilter2.Update();
            //用于将离散的三角面片拼接为连续的等值面
            vtkStripper stripper2 = new vtkStripper();
            stripper2.SetInput(contourFilter2.GetOutput());
            stripper2.Update();
            //Actor：等值线
            vtkDataSetMapper contourMapper2 = new vtkDataSetMapper();
            contourMapper2.SetInput(stripper2.GetOutput());
            vtkActor isolineActor2 = new vtkActor();
            isolineActor2.SetMapper(contourMapper2);

            /*等值线值标注*/
            //构建等值线标注几何数据
            //long linesCount = stripper.GetOutput().GetNumberOfLines(); //等值线的数目
            vtkPoints isolinePoints2 = stripper2.GetOutput().GetPoints(); //等值线的点集
            vtkCellArray isolineCells2 = stripper2.GetOutput().GetLines(); //等值线的线单元数组
            vtkDataArray isolineScalars2 = stripper2.GetOutput().GetPointData().GetScalars(); //等值线的点标量数组
            //等值线标注Label数据
            vtkPolyData labelPolyData2 = new vtkPolyData(); //标注点几何数据
            vtkPoints labelPoints2 = new vtkPoints(); //标注点集
            vtkDoubleArray labelScalars2 = new vtkDoubleArray(); //标注点标量数组
            labelScalars2.SetNumberOfComponents(1);
            //labelScalars.SetName("Isovalues");
            //构建等值线标注Label数据，此处将每一个等值线进行遍历，获取线上的一点作为等值线的标注点
            vtkIdList isoCell2 = new vtkIdList();
            while (isolineCells2.GetNextCell(isoCell2) > 0)
            {
                long idCount2 = isoCell2.GetNumberOfIds(); //获取这条等值线单元的点个数
                //long samplePointIndex = (long)vtkMath.Random(0, idCount); //随机取该等值线上一个点
                long labelPointId2 = isoCell2.GetId(idCount2 / 2); //获取该点的id作为标记点。此方法是根据点数量顺序排序序号Index获取，获取到这个cell中的指定序号Index的点的vtkId（个人认为）

                double[] midPoint2 = isolinePoints2.GetPoint(labelPointId2); //获取该点的坐标

                labelPoints2.InsertNextPoint(midPoint2[0], midPoint2[1], midPoint2[2]); //插入到标注点集
                labelScalars2.InsertNextTuple1(isolineScalars2.GetTuple1(labelPointId2)); //取该点的标量值插入到标注点标量数组
            }
            labelPolyData2.SetPoints(labelPoints2); //设置标注点集
            labelPolyData2.GetPointData().SetScalars(labelScalars2); //设置点数据的标量数组
            labelPolyData2.Update();
            //Mapper
            vtkLabeledDataMapper labelMapper2 = new vtkLabeledDataMapper();
            //labelMapper.SetFieldDataName("Isovalues");
            labelMapper2.SetInput(labelPolyData2);
            labelMapper2.SetLabelModeToLabelScalars();
            labelMapper2.SetLabelFormat("%6.2f");
            labelMapper2.GetLabelTextProperty().SetColor(0, 0, 1);
            //Actor
            vtkActor2D isoLabelActor2 = new vtkActor2D();
            isoLabelActor2.SetMapper(labelMapper2);




            // Create a vtkCamera, and set the camera parameters.
            camera = vtkCamera.New();
            camera.SetClippingRange(1.60187, 20.0842);
            camera.SetFocalPoint(0.21406, 1.5, 0);
            camera.SetPosition(8.3761, 4.94858, 4.12505);
            camera.SetViewUp(0.180325, 0.549245, -0.815974);

            // Create a vtkLight, and set the light parameters.
            light = vtkLight.New();
            light.SetFocalPoint(0.21406, 1.5, 0);
            light.SetPosition(8.3761, 4.94858, 4.12505);

            // Create the Renderers.  Assign them the appropriate viewport coordinates,
            // active camera, and light.
            ren1 = vtkRenderer.New();
            ren1.SetViewport(0, 0, 0.5, 1.0);
            ren1.SetActiveCamera(camera);
            ren1.AddLight(light);

            ren2 = vtkRenderer.New();
            ren2.SetViewport(0.5, 0, 1.0, 1.0);
            ren2.SetActiveCamera(camera);
            ren2.AddLight(light);

            // Create the RenderWindow and RenderWindowInteractor.
            renWin = vtkRenderWindow.New();
            renWin.AddRenderer(ren1);
            renWin.AddRenderer(ren2);
            renWin.SetWindowName("VTK - Cube Axes");
            renWin.SetSize(600, 300);
            iren = vtkRenderWindowInteractor.New();
            iren.SetRenderWindow(renWin);

            ren1.AddActor(surfaceActor);
            ren1.AddActor(isolineActor);
            ren1.AddActor2D(isoLabelActor);
            ren1.AddActor(scalarBarActor);

            ren2.AddActor(surfaceActor3);
            ren2.AddActor(isolineActor2);
            ren2.AddActor2D(isoLabelActor2);

            // Render
            renWin.Render();

            // Set the user method (bound to key 'u')
            iren.Initialize();
            iren.Start();

            // Set up a check for aborting rendering.
            renWin.AbortCheckEvt += new vtkObject.vtkObjectEventHandler(TkCheckAbort);

            //Clean Up  
            deleteAllVTKObjects();




        }

        static vtkBYUReader fohe;
        static vtkPolyDataNormals normals;
        static vtkPolyDataMapper foheMapper;
        static vtkLODActor foheActor;
        static vtkOutlineFilter outline;
        static vtkPolyDataMapper mapOutline;
        static vtkActor outlineActor;
        static vtkCamera camera;
        static vtkLight light;
        static vtkRenderer ren1;
        static vtkRenderer ren2;
        static vtkRenderWindow renWin;
        static vtkRenderWindowInteractor iren;
        static vtkTextProperty tprop;
        static vtkCubeAxesActor2D axes;
        static vtkCubeAxesActor2D axes2;
        static int foo;


        /// <summary>
        /// Callback function for renWin.AbortCheckEvt
        /// </summary>+
        /// 
        public static void TkCheckAbort(vtkObject sender, vtkObjectEventArgs e)
        {
            foo = renWin.GetEventPending();
            if ((foo) != 0)
            {
                renWin.SetAbortRender(1);
            }
        }

        ///<summary>Deletes all static objects created</summary>
        public static void deleteAllVTKObjects()
        {
            //clean up vtk objects
            if (fohe != null) { fohe.Dispose(); }
            if (normals != null) { normals.Dispose(); }
            if (foheMapper != null) { foheMapper.Dispose(); }
            if (foheActor != null) { foheActor.Dispose(); }
            if (outline != null) { outline.Dispose(); }
            if (mapOutline != null) { mapOutline.Dispose(); }
            if (outlineActor != null) { outlineActor.Dispose(); }
            if (camera != null) { camera.Dispose(); }
            if (light != null) { light.Dispose(); }
            if (ren1 != null) { ren1.Dispose(); }
            if (ren2 != null) { ren2.Dispose(); }
            if (renWin != null) { renWin.Dispose(); }
            if (iren != null) { iren.Dispose(); }
            if (tprop != null) { tprop.Dispose(); }
            if (axes != null) { axes.Dispose(); }
            if (axes2 != null) { axes2.Dispose(); }
        }
    }
}
