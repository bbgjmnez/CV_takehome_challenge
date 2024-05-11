#pragma once

#include <bitset>
#include <fstream>
#include <stdexcept>
#include <string>

#include <Eigen/Dense>


class PlyElement {
  std::ofstream &ply_stream_;
  int size_offset_;
  int element_count_ = 0;
  int attr_count_ = 0;
  int attr_index_ = 0;
  bool binary_;
  std::bitset<32> is_list_;
  int list_remaining_ = 0;

  // Constructors are private, accessed by the stream classes via friend
  // declarations.
  PlyElement(std::ofstream &stream, bool binary)
      : ply_stream_(stream), binary_(binary) {}
  ~PlyElement() {
    // CHECK_EQ(attr_index_, 0);
    ply_stream_.seekp(size_offset_);
    ply_stream_ << element_count_ << "\ncomment ";
  }

  friend class PointPlyStream;
  friend class FacePlyStream;

public:
  PlyElement(const PlyElement &) = delete;
  PlyElement(PlyElement &&) = delete;
  PlyElement &operator=(const PlyElement &) = delete;
  PlyElement &operator=(PlyElement &&) = delete;

  void WriteHeader(std::string const &element_name,
                   std::initializer_list<std::string> attributes) {
  ply_stream_ << "element " << element_name << " ";
  size_offset_ = ply_stream_.tellp();
  ply_stream_ << "0\n";
  ply_stream_ << "comment          padding\n";
  int i = 0;
  for (const auto &attr_str : attributes) {
    ply_stream_ << "property " << attr_str << "\n";
    if (attr_str.find("list") == 0) {
      is_list_[i] = true;
    }
    ++i;
  }

  attr_count_ = attributes.size();
                   }

  template <typename T, typename = std::enable_if_t<std::is_scalar<T>::value>>
  PlyElement &operator<<(const T &v) {
    if (attr_index_ == 0 && list_remaining_ == 0) {
      element_count_++;
    }

    ply_stream_.write(reinterpret_cast<const char *>(&v), sizeof(v));
    if (list_remaining_ > 0) {
      list_remaining_--;
    } else if (is_list_.test(attr_index_)) {
      list_remaining_ = v;
    }
    if (list_remaining_ == 0) {
      attr_index_ = (attr_index_ + 1) % attr_count_;
    }
    return *this;
  }

  template <class Derived>
  PlyElement &operator<<(const Eigen::PlainObjectBase<Derived> &mat) {
    if (attr_index_ == 0 && list_remaining_ == 0) {
      element_count_++;
    }
    ply_stream_.write(reinterpret_cast<const char *>(mat.data()),
                      mat.size() * sizeof(typename Derived::Scalar));
    if (list_remaining_ > 0) {
      list_remaining_ -= mat.size();
      assert(list_remaining_ >= 0);
    }
    if (list_remaining_ == 0) {
      attr_index_ = (attr_index_ + mat.size()) % attr_count_;
    }
    return *this;
  }

  template <class Derived>
  PlyElement &operator<<(const Eigen::DenseBase<Derived> &mat) {
    return *this << mat.eval();
  }
};

class PointPlyStream {
  std::ofstream ply_stream_;
  PlyElement point_element_;
  bool binary_;

public:
  PointPlyStream(const std::string &path, bool binary = true)
      : ply_stream_(path, std::ios::binary),
        point_element_(ply_stream_, binary), binary_(binary){};

  bool IsOpen() { return ply_stream_.is_open(); }

  void WriteHeader(std::initializer_list<std::string> attributes,
                   const std::vector<std::string> &comments = {}) {

  ply_stream_ << "ply\n";
  ply_stream_ << "format " << (binary_ ? "binary_little_endian" : "ascii")
              << " 1.0\n";
  for (const auto &comment : comments) {
    ply_stream_ << "comment " << comment << "\n";
  }
  point_element_.WriteHeader("vertex", attributes);
  ply_stream_ << "end_header\n";

}

  template <class... Args> PointPlyStream &operator<<(Args... args) {
    point_element_.operator<<(std::forward<Args...>(args...));
    return *this;
  }
};