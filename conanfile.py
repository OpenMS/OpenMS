from conans import ConanFile, CMake

class OpenMS(ConanFile):
   settings = "os", "compiler", "build_type", "arch"
   requires = [
    "eigen/[>3.0]",
    "boost/[>1.60]",
    "coin-cbc/[>1.0]",
    "coin-cgl/[>0.0]",
    "coin-clp/[>1.0]",
    "coin-osi/[>0.0]",
    "coin-utils/[>1.0]",
    "libsvm/[>300]",
    "hdf5/[>1.0]",
    "sqlite3/[>3.0]",
    "zlib/[>=1.2]",
    "bzip2/[>=1.0]",
    "qt/5.15.5",
    "xerces-c/[>=3.2]"
    ]

   tool_requires = "cmake/[>=3.24]"

   generators = "CMakeDeps"

   default_options = {
      "qt:commercial": False,
      "qt:config": None,
      "qt:device": None,
      "qt:gui": True,
      "qt:multiconfiguration": False,
      "qt:opengl": "desktop",
      "qt:openssl": False,
      "qt:qt3d": False,
      "qt:qtandroidextras": False,
      "qt:qtcharts": False,
      "qt:qtconnectivity": False,
      "qt:qtdatavis3d": False,
      "qt:qtdeclarative": False,
      "qt:qtdoc": False,
      "qt:qtgamepad": False,
      "qt:qtgraphicaleffects": False,
      "qt:qtimageformats": False,
      "qt:qtlocation": False,
      "qt:qtlottie": False,
      "qt:qtmultimedia": False,
      "qt:qtnetworkauth": True,
      "qt:qtpurchasing": False,
      "qt:qtquick3d": False,
      "qt:qtquickcontrols": False,
      "qt:qtquickcontrols2": False,
      "qt:qtquicktimeline": False,
      "qt:qtremoteobjects": False,
      "qt:qtscript": False,
      "qt:qtscxml": False,
      "qt:qtsensors": False,
      "qt:qtserialbus": False,
      "qt:qtserialport": False,
      "qt:qtspeech": False,
      "qt:qtsvg": True,
      "qt:qttools": False,
      "qt:qttranslations": False,
      "qt:qtvirtualkeyboard": False,
      "qt:qtwayland": False,
      "qt:qtwebchannel": False,
      "qt:qtwebengine":False,
      "qt:qtwebglplugin": False,
      "qt:qtwebsockets": False,
      "qt:qtwebview": False,
      "qt:qtxmlpatterns": False,
      "qt:shared": True,
      "qt:widgets": True,
      "qt:with_dbus": False,
      "qt:with_doubleconversion": True,
      "qt:with_freetype": True,
      "qt:with_glib": False,
      "qt:with_harfbuzz": False,
      "qt:with_libjpeg": "libjpeg",
      "qt:with_libpng": True,
      "qt:with_md4c": False,
      "qt:with_mysql": False,
      "qt:with_odbc": False,
      "qt:with_pcre2": False,
      "qt:with_pq": False,
      "qt:with_sqlite3": True,
      "qt:with_vulkan": False,
      "qt:with_zstd": False,
   }


   def configure(self):
      if self.settings.os == "Linux":
         self.default_options["qt:qtx11extras"] = True
         self.default_options["qt:with_fontconfig"] = True
         self.default_options["qt:with_icu"] = False
      elif self.settings.os == "Windows":
         self.default_options["qt:qtwinextras"] = True
      elif self.settings.os == "Macos":
         self.default_options["qt:qtmacextras"] = True

      if self.settings.compiler in ["gcc","clang","intel"] :
         self.settings.compiler.libcxx = "libstd++11"

      if self.settings.compiler == "apple-clang":
         self.settings.compiler.libcxx = "libc++"

   def imports(self):
      self.copy("*", dst="./cmake", src=".", root_package="cmake")
