#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkHexahedron.h>
#include <vtkPolyVertex.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>
#include <vector>
#include <math.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vector>
#include <vtkFieldData.h>

int nEx = 8;
int nEy = 6;
int nEz = 6;
int nSx = 2;
int nSy = 1;
int nSz = 1;



bool asciiOrBinaryVtu = true;
bool display_mesh = false;
bool rotateMeshAccordingToPoinconBenchmark = true;
int decomposition = 3;


using namespace std;

int main(int argc, char *argv[])
{

   // - geometry setting ()
    double lengthUp[] = {4.0, 2.0, 2.0};
    // - decomposition
    int nElSub[] = {nEx,nEy,nEz};
    int nSub[]   = {nSx,nSy,nSz};
    // LOWER BODY
    // - geometry setting ()
    // - decomposition


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                      DEFINITION OF FEM BODY
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
    double shiftUp[3];
    int nElxyz_allUp[3];
    for (int i = 0; i < 3; i++){
        nElxyz_allUp[i] = nElSub[i] * nSub[i];
        shiftUp[i] = 0.5 * lengthUp[i];
    }

    int nEl_FEM = nElxyz_allUp[0] * nElxyz_allUp[1] * nElxyz_allUp[2];
    int nEl_IGA = 1;
    int nEl = nEl_IGA + nEl_FEM;


    double dxyzUp[3];
    for (int i = 0 ; i < 3; i++){
        dxyzUp[i] = lengthUp[i] / nElxyz_allUp[i];
    }

    int nP_FEM = (nElxyz_allUp[0] + 1) * (nElxyz_allUp[1] + 1) * (nElxyz_allUp[2] + 1);





    std::string filename = "iga_fem.vtu";//argv[1];


  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();


// -----------------------------------------------------------------------------------------
    vtkSmartPointer<vtkDoubleArray> _vtkDataArray0 = vtkSmartPointer<vtkDoubleArray>::New();
    _vtkDataArray0->SetName("Knots_0");
    _vtkDataArray0->SetNumberOfComponents(1);
    vector < double > t0( {0.0e+00,  0.e+00,  0.0e+00,  1.0e+00,  1.0e+00,  1.0e+00,
                           0.0e+00,  0.e+00,  0.0e+00,  1.0e+00,  1.0e+00,  1.0e+00,
                           0.0e+00,  0.e+00,  0.0e+00,  1.0e+00,  1.0e+00,  1.0e+00});
    _vtkDataArray0->SetNumberOfTuples(t0.size());
    double tuple[] = {0};
    for (int i = 0; i < t0.size(); i++){
        tuple[0] = t0[i];
        _vtkDataArray0->SetTuple(i, tuple);
    }
    unstructuredGrid->GetFieldData()->AddArray(_vtkDataArray0);
// -----------------------------------------------------------------------------------------
    vtkSmartPointer<vtkDoubleArray> _vtkDataArray1 = vtkSmartPointer<vtkDoubleArray>::New();
    _vtkDataArray1->SetName("Weights_0");
    _vtkDataArray1->SetNumberOfComponents(1);
    vector < double > t1(  {1.0e+00,  1.0e+00,  1.0e+00,
                            1.0e+00,  1.0e+00,  1.0e+00,
                            1.0e+00,  1.0e+00,  1.0e+00,
                            1.0e+00,  1.0e+00,  1.0e+00,
                            1.0e+00,  1.0e+00,  1.0e+00,
                            1.0e+00,  1.0e+00,  1.0e+00,
                            1.0e+00,  1.0e+00,  1.0e+00,
                            1.0e+00,  1.0e+00,  1.0e+00,
                            1.0e+00,  1.0e+00,  1.0e+00});
    _vtkDataArray1->SetNumberOfTuples(t1.size());


    for (int i = 0; i < t1.size(); i++){
        tuple[0] = t1[i];
        _vtkDataArray1->SetTuple(i, tuple);
    }
    unstructuredGrid->GetFieldData()->AddArray(_vtkDataArray1);
// -----------------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray2 = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray2->SetName("PatchData_0");
    _vtkDataArray2->SetNumberOfComponents(1);
    vector < double > t2(  {2,  2,  2,  3,  3,  3,  6,  6,  6});
    _vtkDataArray2->SetNumberOfTuples(t2.size());


    for (int i = 0; i < t2.size(); i++){
        tuple[0] = t2[i];
        _vtkDataArray2->SetTuple(i, tuple);
    }
    unstructuredGrid->GetFieldData()->AddArray(_vtkDataArray2);
// -----------------------------------------------------------------------------------------



    int nP_IGA = t2[3] * t2[4] * t2[5];




// number of IGA nodes
//    int nP_IGA = t2[3] * t2[4] * t2[5];
    vtkSmartPointer<vtkIntArray> _vtkDataArray3 = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray3->SetName("RegionId");
    _vtkDataArray3->SetNumberOfComponents(1);
    vector < int > t3();
    _vtkDataArray3->SetNumberOfTuples(nP_IGA + nP_FEM);

    for (int i = 0; i < nP_IGA; i++){
        tuple[0] = 0;
        _vtkDataArray3->SetTuple(i, tuple);
    }
    for (int i = nP_IGA; i < nP_IGA + nP_FEM; i++){
        tuple[0] = 1;
        _vtkDataArray3->SetTuple(i, tuple);
    }
    unstructuredGrid->GetPointData()->AddArray(_vtkDataArray3);
// -----------------------------------------------------------------------------------------

    double _x,_y,_z;
    int nPoints = 0;

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vector < double > vPointsIGA({
         -4.00000000e+00, -1.00000000e+00, -1.00000000e+00,
         -4.00000000e+00, -1.00000000e+00,  0.00000000e+00,
         -4.00000000e+00, -1.00000000e+00,  1.00000000e+00,
         -4.00000000e+00,  0.00000000e+00, -1.00000000e+00,
         -4.00000000e+00,  0.00000000e+00,  0.00000000e+00,
         -4.00000000e+00,  0.00000000e+00,  1.00000000e+00,
         -4.00000000e+00,  1.00000000e+00, -1.00000000e+00,
         -4.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         -4.00000000e+00,  1.00000000e+00,  1.00000000e+00,
         -2.00000000e+00, -1.00000000e+00, -1.00000000e+00,
         -2.00000000e+00, -1.00000000e+00,  0.00000000e+00,
         -2.00000000e+00, -1.00000000e+00,  1.00000000e+00,
         -2.00000000e+00,  0.00000000e+00, -1.00000000e+00,
         -2.00000000e+00,  0.00000000e+00,  0.00000000e+00,
         -2.00000000e+00,  0.00000000e+00,  1.00000000e+00,
         -2.00000000e+00,  1.00000000e+00, -1.00000000e+00,
         -2.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         -2.00000000e+00,  1.00000000e+00,  1.00000000e+00,
          0.00000000e+00, -1.00000000e+00, -1.00000000e+00,
          0.00000000e+00, -1.00000000e+00,  0.00000000e+00,
          0.00000000e+00, -1.00000000e+00,  1.00000000e+00,
          0.00000000e+00,  0.00000000e+00, -1.00000000e+00,
          0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
          0.00000000e+00,  0.00000000e+00,  1.00000000e+00,
          0.00000000e+00,  1.00000000e+00, -1.00000000e+00,
          0.00000000e+00,  1.00000000e+00,  0.00000000e+00,
          0.00000000e+00,  1.00000000e+00,  1.00000000e+00 });

    for (int i = 0; i <  nP_IGA ; i++){
        _x = vPointsIGA[3 * i + 0];
        _y = vPointsIGA[3 * i + 1];;
        _z = vPointsIGA[3 * i + 2];;
        points->InsertNextPoint(_x, _y, _z);
        nPoints++;
    }

    for (int kk = 0; kk <  nElxyz_allUp[2] + 1; kk++){
        for (int jj = 0; jj <  nElxyz_allUp[1] + 1; jj++){
            for (int ii = 0; ii <  nElxyz_allUp[0] + 1; ii++){
                _x = ii * dxyzUp[0] - shiftUp[0] * 0;
                _y = jj * dxyzUp[1] - shiftUp[1];
                _z = kk * dxyzUp[2] - shiftUp[2];
                points->InsertNextPoint(_x, _y, _z);
                nPoints++;
            }
        }
    }
// Hexahedron




    vtkSmartPointer<vtkIdList> plvx_ids = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> hexa_ids = vtkSmartPointer<vtkIdList>::New();



    vector < int > v_plVert({ 0,   9,  18,
                              3,  12,  21,
                              6,  15,  24,
                              1,  10,  19,
                              4,  13,  22,
                              7,  16,  25,
                              2,  11,  20,
                              5,  14,  23,
                              8,  17,  26});

    for (int i = 0 ; i < v_plVert.size(); i++){
        plvx_ids->InsertId(i,v_plVert[i]);
    }
    unstructuredGrid->Allocate(nEl,1);
    unstructuredGrid->InsertNextCell (VTK_POLY_VERTEX, plvx_ids);


    int ttt[8], _first_set[4];
    std::vector < int > tmp_vec (ttt, ttt + sizeof(ttt)/sizeof(int));
    int nxyUp = (nElxyz_allUp[0] + 1) * (nElxyz_allUp[1] + 1);

    for (int kk = 0; kk <  nElxyz_allUp[2]; kk++){
        for (int jj = 0; jj <  nElxyz_allUp[1]; jj++){
            for (int ii = 0; ii <  nElxyz_allUp[0]; ii++)
            {
                _first_set[0] = ii + 0 + (jj + 0) * (nElxyz_allUp[0] + 1);
                _first_set[1] = ii + 1 + (jj + 0) * (nElxyz_allUp[0] + 1);
                _first_set[2] = ii + 1 + (jj + 1) * (nElxyz_allUp[0] + 1);
                _first_set[3] = ii + 0 + (jj + 1) * (nElxyz_allUp[0] + 1);
                for (int ll = 0; ll < 4; ll++){
                    tmp_vec[ll]     =  nP_IGA + _first_set[ll] + kk * nxyUp;
                    tmp_vec[ll+4]   =  nP_IGA + _first_set[ll] + (kk + 1) * nxyUp;
                }
                for (int ll = 0; ll < 8; ll ++){
                    hexa_ids->InsertId(ll,tmp_vec[ll]);
                }
                unstructuredGrid->InsertNextCell (VTK_HEXAHEDRON, hexa_ids);
            }
        }
    }

//---------------------------

    unstructuredGrid->SetPoints(points);




// PartitionId --------------------------------------------------------------------------------


    vtkSmartPointer<vtkIntArray> _vtkDataArray_PartId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_PartId->SetName("PartitionId");
    _vtkDataArray_PartId->SetNumberOfComponents(1);
    _vtkDataArray_PartId->SetNumberOfTuples(nEl);

    std::cout << "decomposition set-up: " << decomposition << std::endl;


    // FEM
    int cnt = 0;
    int currentIdOfMat = 0;
    tuple[0] = currentIdOfMat;
    _vtkDataArray_PartId->SetTuple(cnt, tuple);
    cnt++;
    currentIdOfMat++;
    for (int KK = 0; KK <  nSub[2]; KK++){
        for (int JJ = 0; JJ <  nSub[1]; JJ++){
            for (int II = 0; II <  nSub[0]; II++){
                for (int kk = KK * nElSub[2] ; kk <  (KK + 1) * nElSub[2]; kk++){
                    for (int jj = JJ * nElSub[1]; jj <  (JJ + 1) * nElSub[1]; jj++){
                        for (int ii = II * nElSub[0]; ii <  (II + 1) * nElSub[0]; ii++){
                            tuple[0] = currentIdOfMat;
                            cnt = 1 + ii + jj * nElxyz_allUp[0] + kk * (nElxyz_allUp[0] * nElxyz_allUp[1]);
                            _vtkDataArray_PartId->SetTuple(cnt, tuple);
                        }
                    }
                }
                currentIdOfMat++;
            }
        }
    }





// MaterialId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_MatId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_MatId->SetName("MaterialId");
    _vtkDataArray_MatId->SetNumberOfComponents(1);
    _vtkDataArray_MatId->SetNumberOfTuples(nEl);
    tuple[0] = 1;
    _vtkDataArray_MatId->SetTuple(0, tuple);
    tuple[0] = 2;
    for (int i = nEl_IGA; i <  nEl_IGA + nEl_FEM; i++){
        _vtkDataArray_MatId->SetTuple(i, tuple);
    }
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_MatId);
// FormulationId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_FormId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_FormId->SetName("FormulationId");
    _vtkDataArray_FormId->SetNumberOfComponents(1);
    _vtkDataArray_FormId->SetNumberOfTuples(nEl);
    tuple[0] = 900;
    _vtkDataArray_FormId->SetTuple(0, tuple);
    tuple[0] = 30;
    for (int i = nEl_IGA; i <  nEl_IGA + nEl_FEM; i++){
        _vtkDataArray_FormId->SetTuple(i, tuple);
    }
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_FormId);
// PieceId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_PieceId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_PieceId->SetName("PieceId");
    _vtkDataArray_PieceId->SetNumberOfComponents(1);
    _vtkDataArray_PieceId->SetNumberOfTuples(nEl);
    tuple[0] = 1;
    _vtkDataArray_PieceId->SetTuple(0, tuple);
    tuple[0] = 2;
    for (int i = nEl_IGA; i <  nEl_IGA + nEl_FEM; i++){
        _vtkDataArray_PieceId->SetTuple(i, tuple);
    }
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_PieceId);
// RegionId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_RegionId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_RegionId->SetName("RegionId");
    _vtkDataArray_RegionId->SetNumberOfComponents(1);
    _vtkDataArray_RegionId->SetNumberOfTuples(nEl);
    tuple[0] = 0;
    _vtkDataArray_RegionId->SetTuple(0, tuple);
    tuple[0] = 1;
    for (int i = nEl_IGA; i <  nEl_IGA + nEl_FEM; i++){
        _vtkDataArray_RegionId->SetTuple(i, tuple);
    }
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_RegionId);

    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_PartId);

    // Write file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());

    if (asciiOrBinaryVtu){
      writer->SetDataModeToAscii();
    }
    else{
      writer->SetDataModeToBinary();
    }


#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(unstructuredGrid);
#else
  writer->SetInputData(unstructuredGrid);
#endif
  writer->Write();

  // Read and display file for verification that it was written correclty
    if (display_mesh){
        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
          vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        reader->SetFileName(filename.c_str());
        reader->Update();

        vtkSmartPointer<vtkDataSetMapper> mapper =
          vtkSmartPointer<vtkDataSetMapper>::New();
        mapper->SetInputConnection(reader->GetOutputPort());

        vtkSmartPointer<vtkActor> actor =
          vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);

        vtkSmartPointer<vtkRenderer> renderer =
          vtkSmartPointer<vtkRenderer>::New();
        vtkSmartPointer<vtkRenderWindow> renderWindow =
          vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
          vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        renderer->AddActor(actor);
        renderer->SetBackground(.3, .6, .3); // Background color green

        renderWindow->Render();
        renderWindowInteractor->Start();
    }

  return EXIT_SUCCESS;
}

