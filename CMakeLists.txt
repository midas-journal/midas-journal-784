INCLUDE(${CMAKE_SOURCE_DIR}/images/IJMacros.txt)

CMAKE_MINIMUM_REQUIRED(VERSION 2.5)

#Change PROJECT_NAME to the name of your project
PROJECT(ImageRegistration2D3D)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)

#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists
#
# The packages can be specified with a version number, for example:
#
# ITK 2.8.1
# ITK 3.2.0
#
# If no version is specified, the most recent release of the package
# will be used.
SET(Required_Packages
ITK
)

#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(Libraries
ITKCommon
ITKIO
ITKBasicFilters
ITKNumerics
)

#the following block of code is an example of how to build an executable in
#cmake.  Unmodified, it will add an executable called "MyExe" to the project.
#MyExe will be built using the files MyClass.h and MyClass.cxx, and it will
#be linked to all the libraries you specified above. 
#You can build more than one executable per project
SET(CurrentExe "DicomSeriesReadNiftiImageWrite")
ADD_EXECUTABLE(${CurrentExe} DicomSeriesReadNiftiImageWrite.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

SET(CurrentExe "ReadResampleWriteNifti")
ADD_EXECUTABLE(${CurrentExe} ReadResampleWriteNifti.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

SET(CurrentExe "getDRRSiddonJacobsRayTracing")
ADD_EXECUTABLE(${CurrentExe} getDRRSiddonJacobsRayTracing.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

SET(CurrentExe "2D3DTwoProjRegistration")
ADD_EXECUTABLE(${CurrentExe} 2D3DTwoProjRegistration.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

SET(CurrentExe "ImageCompare")
ADD_EXECUTABLE(${CurrentExe} ImageCompare.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

#the following line is an example of how to add a test to your project.
#Testname is the title for this particular test.  ExecutableToRun is the
#program which will be running this test.  It can either be a part of this
#project or an external executable.  After that list any args that are needed
#for this test.  Include as many tests as you like.  If your project doesn't have
#any tests you can comment out or delete the following line.
ADD_TEST(CreateProjectionImage getDRRSiddonJacobsRayTracing -rp 0 -iso 255 259 130  -rx -3 -ry 4 -rz 2 -t 5 5 5 -res 0.5 0.5 -size 512 512 -o ${CMAKE_SOURCE_DIR}/images/BoxheadDRR_Disp1_G0.tif ${CMAKE_SOURCE_DIR}/images/BoxheadCTScan1.img)
ADD_TEST(CreateProjectionImage getDRRSiddonJacobsRayTracing -rp 90 -iso 255 259 130  -rx -3 -ry 4 -rz 2 -t 5 5 5 -res 0.5 0.5 -size 512 512 -o ${CMAKE_SOURCE_DIR}/images/BoxheadDRR_Disp1_G90.tif ${CMAKE_SOURCE_DIR}/images/BoxheadCTScan1.img)
ADD_TEST(2D3DTwoProjRegistration 2D3DTwoProjRegistration -res 0.5 0.5 0.5 0.5 -iso 255 259 130 -o ${CMAKE_SOURCE_DIR}/images/BoxheadDRR_Disp1_G0_Reg.tif ${CMAKE_SOURCE_DIR}/images/BoxheadDRR_Disp1_G90_Reg.tif ${CMAKE_SOURCE_DIR}/images/BoxheadDRR_Disp1_G0.tif 0 ${CMAKE_SOURCE_DIR}/images/BoxheadDRR_Disp1_G90.tif 90 ${CMAKE_SOURCE_DIR}/images/BoxheadCTScan1.img)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/BoxheadDRR_Disp1_G0_Reg_Author.tif ${CMAKE_SOURCE_DIR}/images/BoxheadDRR_Disp1_G0_Reg.tif)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/BoxheadDRR_Disp1_G90_Reg_Author.tif ${CMAKE_SOURCE_DIR}/images/BoxheadDRR_Disp1_G90_Reg.tif)


