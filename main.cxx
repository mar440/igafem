#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkHexahedron.h>
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


bool asciiOrBinaryVtu = true;
bool display_mesh = false;
bool rotateMeshAccordingToPoinconBenchmark = true;
int decomposition = 3;
// 1 - problem is supposed to be undecomposed
// 2 - problem is decomposed according bodies (upper and lower)
// 3 - problem is decomposed into 2 subdom:
//      first subdomain contains whole lower and half of upper body,
//      second subdomain is half of upper body
// 4 - decomposition into 2 subdomain:
//      poicon (mobile part) is split into 2 subdomains (0 and 1)
//      and subdom 1 continues to second domain (deformable obstacle)
// x - any other number respects uniform decomp.
//          controled by 'nSubUp' and 'nElSubDwn'





          // UPPER BODY ------------------------------------
          // - geometry setting ()
          double lengthUp[] = {6.0, 6.0, 15.0};
          double radiusUp = 5.0;
          // - decomposition
          int nElSubUp[] = {3,2,3};
          int nSubUp[] = {2,3,4};
          // LOWER BODY
          // - geometry setting ()
          // - decomposition



int main(int argc, char *argv[])
{

    std::string filename = "POINCON_3D.vtu";//argv[1];
    const double PI = 3.14159265359;


    double characterSize = lengthUp[0] + lengthUp[1] + lengthUp[2];


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                      DEFINITION OF UPPER BODY (PART OF BALL)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
    double shiftUp[3];
    int nElxyz_allUp[3];
    for (int i = 0; i < 3; i++){
        nElxyz_allUp[i] = nElSubUp[i] * nSubUp[i];
        shiftUp[i] = 0.5 * lengthUp[i];
    }

    int nAllElementsUp = nElxyz_allUp[0] * nElxyz_allUp[1] * nElxyz_allUp[2];
    double dxyzUp[3];
    for (int i = 0 ; i < 3; i++){
        dxyzUp[i] = lengthUp[i] / nElxyz_allUp[i];
    }
    double _x,_y,_z;
//  Points
    int nPointsUp = 0;
    double weight;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int kk = 0; kk <  nElxyz_allUp[2] + 1; kk++){
        for (int jj = 0; jj <  nElxyz_allUp[1] + 1; jj++){
            for (int ii = 0; ii <  nElxyz_allUp[0] + 1; ii++){
                _x = ii * dxyzUp[0] - shiftUp[0];
                _y = jj * dxyzUp[1] - shiftUp[1];
                _z = kk * dxyzUp[2] - shiftUp[2];
                points->InsertNextPoint(_x, _y, _z);
                nPointsUp++;
            }
        }
    }
// Hexahedron
  vtkSmartPointer<vtkHexahedron> hexa = vtkSmartPointer<vtkHexahedron>::New();
// Cell array - connectivity
  vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();

    int ttt[8], _first_set[4];
    std::vector < int > tmp_vec (ttt, ttt + sizeof(ttt)/sizeof(int));
    int nxyUp = (nElxyz_allUp[0] + 1) * (nElxyz_allUp[1] + 1);

    for (int kk = 0; kk <  nElxyz_allUp[2]; kk++){
        for (int jj = 0; jj <  nElxyz_allUp[1]; jj++){
            for (int ii = 0; ii <  nElxyz_allUp[0]; ii++){
                _first_set[0] = ii + 0 + (jj + 0) * (nElxyz_allUp[0] + 1);
                _first_set[1] = ii + 1 + (jj + 0) * (nElxyz_allUp[0] + 1);
                _first_set[2] = ii + 1 + (jj + 1) * (nElxyz_allUp[0] + 1);
                _first_set[3] = ii + 0 + (jj + 1) * (nElxyz_allUp[0] + 1);
                for (int ll = 0; ll < 4; ll++){
                    tmp_vec[ll]     = _first_set[ll] + kk * nxyUp;
                    tmp_vec[ll+4]   = _first_set[ll] + (kk + 1) * nxyUp;
                }
                for (int ll = 0; ll < 8; ll ++){
                    hexa->GetPointIds()->SetId(ll, tmp_vec[ll]);
                }
                cellArray->InsertNextCell(hexa);
            }
        }
    }
   //
    int nPoints = nPointsUp;

//---------------------------

  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  unstructuredGrid->SetPoints(points);
  unstructuredGrid->SetCells(VTK_HEXAHEDRON, cellArray);

// PartitionId --------------------------------------------------------------------------------


//    vtkDataArray *_vtkDataArray_PartId;
//    _vtkDataArray_PartId = vtkDoubleArray::New();
    vtkSmartPointer<vtkIntArray> _vtkDataArray_PartId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_PartId->SetName("PartitionId");
    _vtkDataArray_PartId->SetNumberOfComponents(1);
    _vtkDataArray_PartId->SetNumberOfTuples(nAllElementsUp);





    double tuple[] = {0};


    std::cout << "decomposition set-up: " << decomposition << std::endl;


    // Upper body
    int cnt = 0;
    int currentIdOfMat = 0;
    for (int KK = 0; KK <  nSubUp[2]; KK++){
        for (int JJ = 0; JJ <  nSubUp[1]; JJ++){
            for (int II = 0; II <  nSubUp[0]; II++){
                for (int kk = KK * nElSubUp[2] ; kk <  (KK + 1) * nElSubUp[2]; kk++){
                    for (int jj = JJ * nElSubUp[1]; jj <  (JJ + 1) * nElSubUp[1]; jj++){
                        for (int ii = II * nElSubUp[0]; ii <  (II + 1) * nElSubUp[0]; ii++){
                            tuple[0] = currentIdOfMat;
                            cnt = ii + jj * nElxyz_allUp[0] + kk * (nElxyz_allUp[0] * nElxyz_allUp[1]);
                            _vtkDataArray_PartId->SetTuple(cnt, tuple);
                        }
                    }
                }
                currentIdOfMat++;
            }
        }
    }





// MaterialId --------------------------------------------------------------------------------
//
//    vtkDataArray *_vtkDataArray_MatId;
//    _vtkDataArray_MatId = vtkDoubleArray::New();
    vtkSmartPointer<vtkIntArray> _vtkDataArray_MatId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_MatId->SetName("MaterialId");
    _vtkDataArray_MatId->SetNumberOfComponents(1);
    _vtkDataArray_MatId->SetNumberOfTuples(nAllElementsUp);


    tuple[0] = 0;
    for (int i = 0; i <  nAllElementsUp; i++){
        _vtkDataArray_MatId->SetTuple(i, tuple);
    }
    tuple[0] = 1;
    for (int i = nAllElementsUp; i <  nAllElementsUp; i++){
        _vtkDataArray_MatId->SetTuple(i, tuple);
    }
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_MatId);

//
// FormulationId --------------------------------------------------------------------------------
//
//    vtkDataArray *_vtkDataArray_FormId;
//    _vtkDataArray_FormId = vtkDoubleArray::New();
    vtkSmartPointer<vtkIntArray> _vtkDataArray_FormId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_FormId->SetName("FormulationId");
    _vtkDataArray_FormId->SetNumberOfComponents(1);
    _vtkDataArray_FormId->SetNumberOfTuples(nAllElementsUp);


    tuple[0] = 0;
    for (int i = 0; i <  nAllElementsUp; i++){
        _vtkDataArray_FormId->SetTuple(i, tuple);
    }
    tuple[0] = 1;
    for (int i = nAllElementsUp; i <  nAllElementsUp; i++){
        _vtkDataArray_FormId->SetTuple(i, tuple);
    }
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_FormId);
//
// PieceId --------------------------------------------------------------------------------
//
//    vtkDataArray *_vtkDataArray_PieceId;
//    _vtkDataArray_PieceId = vtkDoubleArray::New();
    vtkSmartPointer<vtkIntArray> _vtkDataArray_PieceId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_PieceId->SetName("PieceId");
    _vtkDataArray_PieceId->SetNumberOfComponents(1);
    _vtkDataArray_PieceId->SetNumberOfTuples(nAllElementsUp);


    tuple[0] = 0;
    for (int i = 0; i <  nAllElementsUp; i++){
        _vtkDataArray_PieceId->SetTuple(i, tuple);
    }
    tuple[0] = 1;
    for (int i = nAllElementsUp; i <  nAllElementsUp; i++){
        _vtkDataArray_PieceId->SetTuple(i, tuple);
    }
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_PieceId);

// ---------------


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

