/*
* itkio - a cross-platform C++ library for simplified IO of image data using the ITK library.
*
* Copyright 2016 Ben Glocker <b.glocker@imperial.ac.uk>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#include "itkio.h"

#include "boost/filesystem.hpp"

#include "miaImageIterator.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkImportImageFilter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

namespace itkio
{
  template <typename PixelType>
  void save(const mia::Image& image, const std::string& filename)
  {    
    auto ext = boost::filesystem::extension(filename);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if (image.sizeZ() == 1 && (ext == ".png" || ext == ".jpg" || ext == ".jpeg"))
    {
      const unsigned int Dimension = 2;
      typedef typename itk::Image<PixelType, Dimension> ImageType;
      typedef typename itk::ImportImageFilter<PixelType, Dimension> ImportFilterType;

      typename ImportFilterType::Pointer importer = ImportFilterType::New();

      typename ImageType::SizeType size;
      size[0] = image.sizeX();
      size[1] = image.sizeY();

      typename ImageType::IndexType start;
      start.Fill(0);

      typename ImageType::RegionType region;
      region.SetSize(size);
      region.SetIndex(start);
      importer->SetRegion(region);

      typename ImageType::SpacingType spacing;
      spacing[0] = image.spacing()[0];
      spacing[1] = image.spacing()[1];
      importer->SetSpacing(spacing);

      typename ImageType::PointType origin;
      origin[0] = image.origin()[0];
      origin[1] = image.origin()[1];
      importer->SetOrigin(origin);

      typename ImageType::DirectionType direction;
      for (int r = 0; r < Dimension; r++)
      {
        for (int c = 0; c < Dimension; c++)
        {
          direction(r, c) = image.imageToWorld()(r, c);
        }
      }
      importer->SetDirection(direction);

      long counter = 0;
      std::vector<PixelType> dataConvert(image.size());
      for (auto& p : const_cast<mia::Image&>(image))
      {
        dataConvert[counter++] = static_cast<PixelType>(p);
      }

      importer->SetImportPointer(&dataConvert[0], image.size(), false);
      importer->Update();

      typedef typename  itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName(filename);
      writer->SetInput(importer->GetOutput());
      writer->Update();
    }
    else
    {
      const unsigned int Dimension = 3;
      typedef typename itk::Image<PixelType, Dimension> ImageType;
      typedef typename itk::ImportImageFilter<PixelType, Dimension> ImportFilterType;

      typename ImportFilterType::Pointer importer = ImportFilterType::New();

      typename ImageType::SizeType size;
      size[0] = image.sizeX();
      size[1] = image.sizeY();
      size[2] = image.sizeZ();

      typename ImageType::IndexType start;
      start.Fill(0);

      typename ImageType::RegionType region;
      region.SetSize(size);
      region.SetIndex(start);
      importer->SetRegion(region);

      typename ImageType::SpacingType spacing;
      spacing[0] = image.spacing()[0];
      spacing[1] = image.spacing()[1];
      spacing[2] = image.spacing()[2];
      importer->SetSpacing(spacing);

      typename ImageType::PointType origin;
      origin[0] = image.origin()[0];
      origin[1] = image.origin()[1];
      origin[2] = image.origin()[2];
      importer->SetOrigin(origin);

      typename ImageType::DirectionType direction;
      for (int r = 0; r < Dimension; r++)
      {
        for (int c = 0; c < Dimension; c++)
        {
          direction(r, c) = image.imageToWorld()(r, c);
        }
      }
      importer->SetDirection(direction);

      long counter = 0;
      std::vector<PixelType> dataConvert(image.size());
      for (auto& p : const_cast<mia::Image&>(image))
      {
        dataConvert[counter++] = static_cast<PixelType>(p);
      }

      importer->SetImportPointer(&dataConvert[0], image.size(), false);
      importer->Update();

      if (ext == ".dcm" || ext == ".ima" || ext == "")
      {
        auto pos_ext = filename.rfind(ext);
        auto fname = filename;
        fname.erase(pos_ext, ext.length());
        std::vector<std::string> filenames;
        for (int i = 0; i < image.sizeZ(); i++)
        {

          std::stringstream filename_slice;
          filename_slice << fname << "_" << (i + 1000000) << ext;
          filenames.push_back(filename_slice.str());
        }

        typedef itk::GDCMImageIO ImageIOType;
        ImageIOType::Pointer dicomIO = ImageIOType::New();
        //for (int d = 0; d < 3; d++)
        //{
        //  dicomIO->SetSpacing(d, spacing[d]);
        //  //dicomIO->SetDirection(d, );
        //  dicomIO->SetOrigin(d, origin[d]);
        //}

        typedef typename itk::Image< PixelType, 3> ImageType2D;
        typedef typename  itk::ImageSeriesWriter<ImageType, ImageType2D> WriterType;
        typename WriterType::Pointer writer = WriterType::New();

        writer->SetImageIO(dicomIO);
        writer->SetFileNames(filenames);
        writer->SetInput(importer->GetOutput());
        writer->Update();

        //for (int i = 0; i < image.sizeZ(); i++)
        //{
        //  typedef itk::GDCMImageIO ImageIOType;
        //  ImageIOType::Pointer dicomIO2D = ImageIOType::New();
        //  typedef itk::ImageFileReader< ImageType2D > ReaderType;
        //  ReaderType::Pointer reader2D = ReaderType::New();
        //  reader2D->SetImageIO(dicomIO2D);
        //  reader2D->SetFileName(filenames[i]);
        //  reader2D->Update();

        //  ImageType2D::Pointer inputImage = reader2D->GetOutput();
        //  typedef itk::MetaDataDictionary DictionaryType;
        //  DictionaryType & dictionary = inputImage->GetMetaDataDictionary();
        //  //itk::EncapsulateMetaData<std::string>(dictionary, "0018|0050", std::string("5.0"));
        //  std::stringstream imagepos;
        //  imagepos << std::to_string(origin[0]) << "\\" << std::to_string(origin[1]) << "\\" << std::to_string(origin[2]);
        //  itk::EncapsulateMetaData<std::string>(dictionary, "0020|0032", imagepos.str());

        //  typedef itk::ImageFileWriter< ImageType2D > WriterType2D;
        //  WriterType2D::Pointer writer2D = WriterType2D::New();
        //  writer2D->SetInput(reader2D->GetOutput());
        //  writer2D->SetFileName(filenames[i]);
        //  writer2D->SetImageIO(dicomIO2D);
        //  writer2D->Update();
        //}
      }
      else
      {
        typedef typename  itk::ImageFileWriter<ImageType> WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetFileName(filename);
        writer->SetInput(importer->GetOutput());
        writer->Update();
      }
    }    
  }

  void save(const mia::Image& image, const std::string& filename)
  {
    switch(image.dataType())
    {
    case mia::ImageDataType::BYTE:
      itkio::save<unsigned char>(image, filename);
      break;
    case mia::ImageDataType::SHORT:
      itkio::save<short>(image, filename);
      break;
    case mia::ImageDataType::USHORT:
      itkio::save<unsigned short>(image, filename);
      break;
    case mia::ImageDataType::INT:
      itkio::save<int>(image, filename);
      break;
    case mia::ImageDataType::UINT:
      itkio::save<unsigned char>(image, filename);
      break;
    case mia::ImageDataType::DOUBLE:
      itkio::save<double>(image, filename);
      break;
    default:
      itkio::save<float>(image, filename);
      break;
    }
  }

  mia::Image load(const std::string& filename)
  {
    auto ext = boost::filesystem::extension(filename);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
	if (ext != ".gz" && ext != ".nii" && ext != ".nrrd" && ext != ".mhd" && ext != ".mha" && ext != ".hdr" && ext != ".png" && ext != ".jpg" && ext != ".jpeg")
    {
      auto dicom_image = load_dicom(filename);
      if (dicom_image.size() > 0) return dicom_image;
    }

    typedef float PixelType;
    const unsigned int Dimension = 3;
    typedef itk::Image< PixelType, Dimension > ImageType;
    typedef itk::ImageFileReader< ImageType > ReaderType;

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(filename);
    reader->UpdateOutputInformation();

    ImageType::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();

    itk::ImageIOBase::IOComponentType componentType = reader->GetImageIO()->GetComponentType();
    ComponentType compType = static_cast<ComponentType>(componentType);

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(filename.c_str(), itk::ImageIOFactory::ReadMode);
    imageIO->SetFileName(filename.c_str());
    imageIO->ReadImageInformation();

    int sizeX = imageIO->GetDimensions(0);
    int sizeY = imageIO->GetDimensions(1);
    int sizeZ = imageIO->GetNumberOfDimensions() > 2 ? imageIO->GetDimensions(2) : 1;
    double spacingX = imageIO->GetSpacing(0);
    double spacingY = imageIO->GetSpacing(1);
    double spacingZ = imageIO->GetNumberOfDimensions() > 2 ? imageIO->GetSpacing(2) : 1;
    double originX = imageIO->GetOrigin(0);
    double originY = imageIO->GetOrigin(1);
    double originZ = imageIO->GetNumberOfDimensions() > 2 ? imageIO->GetOrigin(2) : 0;

    mia::Image image(sizeX, sizeY, sizeZ);
    image.spacing(Eigen::Vector3d(spacingX, spacingY, spacingZ));
    image.origin(Eigen::Vector3d(originX, originY, originZ));

    Eigen::Matrix3d imageToWorld = Eigen::Matrix3d::Identity();
    if (sizeZ > 1)
    {
      for (int r = 0; r < Dimension; r++)
      {
        for (int c = 0; c < Dimension; c++)
        {
          imageToWorld(r,c) = imageIO->GetDirection(c)[r];
        }
      }
    }
    image.imageToWorld(imageToWorld);

    typedef itk::ImportImageFilter<PixelType, Dimension> ImportFilterType;

    ImportFilterType::Pointer importer = ImportFilterType::New();

    ImageType::IndexType start;
    start.Fill(0);

    region.SetIndex(start);

    importer->SetRegion(region);

    importer->SetImportPointer(image.data(), image.size(), false);
    importer->Update();

    ImageType::Pointer itkImage = importer->GetOutput();

    reader->GraftOutput(itkImage);
    reader->Update();

    switch(compType)
    {
    case UCHAR:
      image.dataType(mia::ImageDataType::BYTE);
      break;
    case CHAR:
      image.dataType(mia::ImageDataType::BYTE);
      break;
    case SHORT:
      image.dataType(mia::ImageDataType::SHORT);
      break;
    case USHORT:
      image.dataType(mia::ImageDataType::USHORT);
      break;
    case INT:
      image.dataType(mia::ImageDataType::INT);
      break;
    case UINT:
      image.dataType(mia::ImageDataType::UINT);
      break;
    case FLOAT:
      image.dataType(mia::ImageDataType::FLOAT);
      break;
    case DOUBLE:
      image.dataType(mia::ImageDataType::DOUBLE);
      break;
    default:
      image.dataType(mia::ImageDataType::FLOAT);
      break;
    }

    return image;
  }

  mia::Image load_dicom(const std::string& filename)
  {
    boost::filesystem::path p(filename);
    auto directory = p.parent_path();

    const unsigned int Dimension = 3;
    typedef float PixelType;
    typedef itk::Image<PixelType, Dimension> ImageType;
    typedef itk::ImageSeriesReader<ImageType> ReaderType;
    typedef itk::GDCMImageIO ImageIOType;

    itk::GDCMSeriesFileNames::Pointer generator = itk::GDCMSeriesFileNames::New();
    generator->SetUseSeriesDetails(true);

    generator->AddSeriesRestriction("0008|0021"); //series data
    generator->AddSeriesRestriction("0008|0016"); //sop class uid
    generator->AddSeriesRestriction("0020|000E"); //series instance uid
    generator->AddSeriesRestriction("0020|0037"); //image orientation patient
    generator->AddSeriesRestriction("0028|0030"); //pixel spacing
    generator->AddSeriesRestriction("0018|1164"); //imager pixel spacing
    generator->AddSeriesRestriction("0018|0050"); //slice thickness
    generator->AddSeriesRestriction("0028|0010"); //rows
    generator->AddSeriesRestriction("0028|0011"); //columns
    generator->AddSeriesRestriction("0018|1120"); //gantry detector tilt
    generator->AddSeriesRestriction("0008|0060"); //modality
    generator->AddSeriesRestriction("0028|0008"); //number of frames
    generator->AddSeriesRestriction("0018|1140"); //rotation direction
    generator->AddSeriesRestriction("0018|6024"); //physical units x direction
    generator->AddSeriesRestriction("0018|6026"); //physical units y direction
    generator->AddSeriesRestriction("0054|0500"); //slice progression direction
    generator->AddSeriesRestriction("0072|0604"); //sorting direction

    auto filename_lower = p.filename().string();
    std::transform(filename_lower.begin(), filename_lower.end(), filename_lower.begin(), ::tolower);
    generator->SetRecursive(false);
    generator->SetDirectory(directory.string());
    std::vector<std::string> uids = generator->GetSeriesUIDs();

    if (uids.size() > 0)
    {
      std::vector<std::string> filenames = generator->GetFileNames(uids.front());
      for (int id = 0; id < uids.size(); id++)
      {
        filenames.clear();
        filenames = generator->GetFileNames(uids[id]);

        bool stopSearch = false;
        for (int f = 0; f < filenames.size(); f++)
        {
          auto p_search = boost::filesystem::path(filenames[f]);
          auto filename_search_lower = p_search.filename().string();
          std::transform(filename_search_lower.begin(), filename_search_lower.end(), filename_search_lower.begin(), ::tolower);
          stopSearch = filename_lower == filename_search_lower;
          if (stopSearch) break;
        }
        if (stopSearch) break;
      }

      ReaderType::Pointer reader = ReaderType::New();
      ImageIOType::Pointer dicomIO = ImageIOType::New();

      reader->SetImageIO(dicomIO);
      reader->SetFileNames(filenames);

      reader->UpdateOutputInformation();

      ImageType::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();

      itk::ImageIOBase::IOComponentType componentType = reader->GetImageIO()->GetComponentType();
      ComponentType compType = static_cast<ComponentType>(componentType);

      ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
      ImageType::PointType origin = reader->GetOutput()->GetOrigin();
      ImageType::RegionType::SizeType size = region.GetSize();
      ImageType::DirectionType direction = reader->GetOutput()->GetDirection();

      mia::Image image(size[0], size[1], size[2]);
      image.spacing(Eigen::Vector3d(spacing[0], spacing[1], spacing[2]));
      image.origin(Eigen::Vector3d(origin[0], origin[1], origin[2]));

      Eigen::Matrix3d imageToWorld = Eigen::Matrix3d::Identity();
      if (size[2] > 1)
      {
        for (int r = 0; r < Dimension; r++)
        {
          for (int c = 0; c < Dimension; c++)
          {
            imageToWorld(r, c) = direction(r, c);
          }
        }
      }
      image.imageToWorld(imageToWorld);

      typedef itk::ImportImageFilter< PixelType, Dimension > ImportFilterType;

      ImportFilterType::Pointer importer = ImportFilterType::New();

      ImageType::IndexType start;
      start.Fill(0);

      region.SetIndex(start);

      importer->SetRegion(region);

      importer->SetImportPointer(image.data(), image.size(), false);
      importer->Update();

      ImageType::Pointer itkImage = importer->GetOutput();

      reader->GraftOutput(itkImage);
      reader->Update();

      switch (compType)
      {
      case UCHAR:
        image.dataType(mia::ImageDataType::BYTE);
        break;
      case CHAR:
        image.dataType(mia::ImageDataType::BYTE);
        break;
      case SHORT:
        image.dataType(mia::ImageDataType::SHORT);
        break;
      case USHORT:
        image.dataType(mia::ImageDataType::USHORT);
        break;
      case INT:
        image.dataType(mia::ImageDataType::INT);
        break;
      case UINT:
        image.dataType(mia::ImageDataType::UINT);
        break;
      case FLOAT:
        image.dataType(mia::ImageDataType::FLOAT);
        break;
      case DOUBLE:
        image.dataType(mia::ImageDataType::DOUBLE);
        break;
      default:
        image.dataType(mia::ImageDataType::FLOAT);
        break;
      }

      return image;
    }

    return mia::Image();
  }

  std::vector<std::string> get_uids(const std::string& filename)
  {
    boost::filesystem::path p(filename);
    auto directory = p.parent_path();

    const unsigned int Dimension = 3;

    typedef float PixelType;
    typedef itk::Image<PixelType, Dimension> ImageType;
    typedef itk::ImageSeriesReader<ImageType> ReaderType;
    typedef itk::GDCMImageIO ImageIOType;

    itk::GDCMSeriesFileNames::Pointer generator = itk::GDCMSeriesFileNames::New();
    generator->SetUseSeriesDetails(true);

    generator->AddSeriesRestriction("0008|0021"); //series data
    generator->AddSeriesRestriction("0008|0016"); //sop class uid
    generator->AddSeriesRestriction("0020|000E"); //series instance uid
    generator->AddSeriesRestriction("0020|0037"); //image orientation patient
    generator->AddSeriesRestriction("0028|0030"); //pixel spacing
    generator->AddSeriesRestriction("0018|1164"); //imager pixel spacing
    generator->AddSeriesRestriction("0018|0050"); //slice thickness
    generator->AddSeriesRestriction("0028|0010"); //rows
    generator->AddSeriesRestriction("0028|0011"); //columns
    generator->AddSeriesRestriction("0018|1120"); //gantry detector tilt
    generator->AddSeriesRestriction("0008|0060"); //modality
    generator->AddSeriesRestriction("0028|0008"); //number of frames
    generator->AddSeriesRestriction("0018|1140"); //rotation direction
    generator->AddSeriesRestriction("0018|6024"); //physical units x direction
    generator->AddSeriesRestriction("0018|6026"); //physical units y direction
    generator->AddSeriesRestriction("0054|0500"); //slice progression direction
    generator->AddSeriesRestriction("0072|0604"); //sorting direction

    generator->SetRecursive(false);
    generator->SetDirectory(directory.string());
    std::vector<std::string> uids = generator->GetSeriesUIDs();

    std::vector<std::string> uidArray(uids.size());
    for (int id = 0; id < uids.size(); id++)
    {
      uidArray[id] = uids[id].c_str();
    }

    return uidArray;
  }

  itk::Image<float, 3>::Pointer miaToItk(const mia::Image& image)
  {
    const unsigned int Dimension = 3;
    typedef itk::Image<float, Dimension> ImageType;
    typedef itk::ImportImageFilter<float, Dimension> ImportFilterType;

    ImportFilterType::Pointer importer = ImportFilterType::New();

    ImageType::SizeType size;
    size[0] = image.sizeX();
    size[1] = image.sizeY();
    size[2] = image.sizeZ();

    ImageType::IndexType start;
    start.Fill(0);

    ImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);
    importer->SetRegion(region);

    ImageType::SpacingType spacing;
    spacing[0] = image.spacing()[0];
    spacing[1] = image.spacing()[1];
    spacing[2] = image.spacing()[2];
    importer->SetSpacing(spacing);

    ImageType::PointType origin;
    origin[0] = image.origin()[0];
    origin[1] = image.origin()[1];
    origin[2] = image.origin()[2];
    importer->SetOrigin(origin);

    ImageType::DirectionType direction;
    for (int r = 0; r < Dimension; r++)
    {
      for (int c = 0; c < Dimension; c++)
      {
        direction(r, c) = image.imageToWorld()(r, c);
      }
    }
    importer->SetDirection(direction);

    importer->SetImportPointer(const_cast<float*>(image.data()), image.size(), false);
    importer->Update();

    return importer->GetOutput();
  }
}
