ENV_PREFIX = $(CONDA_PREFIX)

CXX = $(ENV_PREFIX)/bin/clang++

CXXFLAGS += \
  -stdlib=libc++ \
  -isystem   '$(ENV_PREFIX)/include/c++/v1' \
  -isystem   '$(ENV_PREFIX)/include' \
  -L         '$(ENV_PREFIX)/lib' \
  -Wl,-rpath,'$(ENV_PREFIX)/lib' \
  -Wno-unused-command-line-argument

OPTFLAGS += \
  -march=x86-64-v3 \
  -mtune=znver2 \
  -flto \
  -fuse-ld=lld

DBGFLAGS += \
  -DNDEBUG_
