/*
* Copyright (c) 2014, Autonomous Systems Lab
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
* * Redistributions in binary form must reproduce the above copyright
* notice, this list of conditions and the following disclaimer in the
* documentation and/or other materials provided with the distribution.
* * Neither the name of the Autonomous Systems Lab, ETH Zurich nor the
* names of its contributors may be used to endorse or promote products
* derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/

#ifndef BAG_UTILS_HPP_
#define BAG_UTILS_HPP_

#include <sstream>

#include <sensor_msgs/Image.h>
#include <sensor_msgs/CompressedImage.h>
#include <sensor_msgs/image_encodings.h>
#include <cv_bridge/cv_bridge.h>

#include <opencv2/imgproc/imgproc.hpp>

namespace rovio{

namespace enc = sensor_msgs::image_encodings;

sensor_msgs::ImageConstPtr decompressImageMsg(const sensor_msgs::CompressedImageConstPtr& message) {
  // Decompresses a CompressedImage to an Image message, according to:
  // https://github.com/ros-perception/image_transport_plugins/blob/indigo-devel/compressed_image_transport/src/compressed_subscriber.cpp

  namespace enc = sensor_msgs::image_encodings;

  cv_bridge::CvImagePtr cv_ptr(new cv_bridge::CvImage);

  // Copy message header
  cv_ptr->header = message->header;

  // Decode color/mono image
  try
  {
    cv_ptr->image = cv::imdecode(cv::Mat(message->data), cv::IMREAD_UNCHANGED);

    // Assign image encoding string
    const size_t split_pos = message->format.find(';');
    if (split_pos==std::string::npos)
    {
      // Older version of compressed_image_transport does not signal image format
      switch (cv_ptr->image.channels())
      {
        case 1:
          cv_ptr->encoding = enc::MONO8;
          break;
        case 3:
          cv_ptr->encoding = enc::BGR8;
          break;
        default:
          ROS_ERROR("Unsupported number of channels: %i", cv_ptr->image.channels());
          break;
      }
    } else
    {
      std::string image_encoding = message->format.substr(0, split_pos);

      cv_ptr->encoding = image_encoding;

      if ( enc::isColor(image_encoding))
      {
        std::string compressed_encoding = message->format.substr(split_pos);
        bool compressed_bgr_image = (compressed_encoding.find("compressed bgr") != std::string::npos);

        // Revert color transformation
        if (compressed_bgr_image)
        {
          // if necessary convert colors from bgr to rgb
          if ((image_encoding == enc::RGB8) || (image_encoding == enc::RGB16))
            cv::cvtColor(cv_ptr->image, cv_ptr->image, CV_BGR2RGB);

          if ((image_encoding == enc::RGBA8) || (image_encoding == enc::RGBA16))
            cv::cvtColor(cv_ptr->image, cv_ptr->image, CV_BGR2RGBA);

          if ((image_encoding == enc::BGRA8) || (image_encoding == enc::BGRA16))
            cv::cvtColor(cv_ptr->image, cv_ptr->image, CV_BGR2BGRA);
        } else
        {
          // if necessary convert colors from rgb to bgr
          if ((image_encoding == enc::BGR8) || (image_encoding == enc::BGR16))
            cv::cvtColor(cv_ptr->image, cv_ptr->image, CV_RGB2BGR);

          if ((image_encoding == enc::BGRA8) || (image_encoding == enc::BGRA16))
            cv::cvtColor(cv_ptr->image, cv_ptr->image, CV_RGB2BGRA);

          if ((image_encoding == enc::RGBA8) || (image_encoding == enc::RGBA16))
            cv::cvtColor(cv_ptr->image, cv_ptr->image, CV_RGB2RGBA);
        }
      }
    }
  }
  catch (cv::Exception& e)
  {
    ROS_ERROR("%s", e.what());
  }

  size_t rows = cv_ptr->image.rows;
  size_t cols = cv_ptr->image.cols;

  if ((rows > 0) && (cols > 0)) {
    return cv_ptr->toImageMsg();
  } else {
    return nullptr;
  }
}

sensor_msgs::ImageConstPtr bagIteratorToImg(const rosbag::View::iterator& it) {
  // Instantiates an Image from the bag iterator, decompressing the
  // payload if needed.
  if (it->getDataType() == "sensor_msgs/CompressedImage") {
    sensor_msgs::CompressedImageConstPtr msg = it->instantiate<sensor_msgs::CompressedImage>();
    return decompressImageMsg(msg);
  } else {
    return it->instantiate<sensor_msgs::Image>();
  }
}

sensor_msgs::CompressedImage compressImageMsg(const sensor_msgs::Image &message) {
  // compresses an Image message to CompressedImage

  sensor_msgs::CompressedImage compressed;
  compressed.header = message.header;
  compressed.format = message.encoding;
  compressed.format += "; jpeg compressed";

  // Compression settings
  std::vector<int> params = {
      CV_IMWRITE_JPEG_QUALITY,
      80
  };
  params.resize(3, 0);

  int bitDepth = enc::bitDepth(message.encoding);
  int numChannels = enc::numChannels(message.encoding);

  // Target image format
  std::stringstream targetFormat;
  if (enc::isColor(message.encoding)) {
    // convert color images to RGB domain
    targetFormat << "rgb" << bitDepth;
  }

  cv_bridge::CvImagePtr cv_ptr;
  try {
    cv_ptr = cv_bridge::toCvCopy(message, targetFormat.str());

    // Compress image
    if (cv::imencode(".jpg", cv_ptr->image, compressed.data, params)) {
      float cRatio = (float) (cv_ptr->image.rows * cv_ptr->image.cols * cv_ptr->image.elemSize())
          / (float) compressed.data.size();
      ROS_DEBUG("Compressed Image Transport - Codec: jpg, Compression: 1:%.2f (%lu bytes)",
                cRatio,
                compressed.data.size());
      return compressed;
    } else {
      ROS_ERROR("cv::imencode (jpeg) failed on input image");
    }
  } catch (cv_bridge::Exception &e) {
    ROS_ERROR("%s", e.what());
  } catch (cv::Exception &e) {
    ROS_ERROR("%s", e.what());
  }
}

}


#endif /* BAG_UTILS_HPP_ */
