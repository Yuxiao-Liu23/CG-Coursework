PROJECT_NAME := renderer-50
BUILD_DIR := build

SOURCE_FILE := $(PROJECT_NAME).cpp
OBJECT_FILE := $(BUILD_DIR)/$(PROJECT_NAME).o
EXECUTABLE := $(BUILD_DIR)/$(PROJECT_NAME)

# sdw 库配置
SDW_DIR := ./sdw/
GLM_DIR := ./libs/glm-0.9.7.2/
SDW_SOURCE_FILES := $(wildcard $(SDW_DIR)*.cpp)
SDW_OBJECT_FILES := $(patsubst $(SDW_DIR)%.cpp, $(BUILD_DIR)/%.o, $(SDW_SOURCE_FILES))

# 编译器设置
COMPILER := clang++
COMPILER_OPTIONS := -c -pipe -Wall -std=c++17
DEBUG_OPTIONS := -ggdb -g3
FUSSY_OPTIONS := -Werror -pedantic
SANITIZER_OPTIONS := -O1 -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer
SPEEDY_OPTIONS := -Ofast -funsafe-math-optimizations -march=native
LINKER_OPTIONS :=

# 编译标志
SDW_COMPILER_FLAGS := -I.
GLM_COMPILER_FLAGS := -I$(GLM_DIR)
SDL_COMPILER_FLAGS := $(shell sdl2-config --cflags 2>/dev/null)
SDL_LINKER_FLAGS := $(shell sdl2-config --libs 2>/dev/null)
SDW_LINKER_FLAGS := $(SDW_OBJECT_FILES)

# SDL2 fallback
ifeq ($(SDL_COMPILER_FLAGS),)
SDL_COMPILER_FLAGS := -I/usr/include/SDL2 -D_REENTRANT
endif
ifeq ($(SDL_LINKER_FLAGS),)
SDL_LINKER_FLAGS := -lSDL2
endif

default: debug

# Debug 模式（推荐用于开发）
debug: $(SDW_OBJECT_FILES)
	@mkdir -p $(BUILD_DIR)
	$(COMPILER) $(COMPILER_OPTIONS) $(DEBUG_OPTIONS) -o $(OBJECT_FILE) $(SOURCE_FILE) $(SDL_COMPILER_FLAGS) $(SDW_COMPILER_FLAGS) $(GLM_COMPILER_FLAGS)
	$(COMPILER) $(LINKER_OPTIONS) $(DEBUG_OPTIONS) -o $(EXECUTABLE) $(OBJECT_FILE) $(SDW_LINKER_FLAGS) $(SDL_LINKER_FLAGS)
	@echo "编译完成！运行: ./$(EXECUTABLE)"
	./$(EXECUTABLE)

# 诊断模式（用于查找内存错误）
diagnostic: $(SDW_OBJECT_FILES)
	@mkdir -p $(BUILD_DIR)
	$(COMPILER) $(COMPILER_OPTIONS) $(FUSSY_OPTIONS) $(SANITIZER_OPTIONS) -o $(OBJECT_FILE) $(SOURCE_FILE) $(SDL_COMPILER_FLAGS) $(SDW_COMPILER_FLAGS) $(GLM_COMPILER_FLAGS)
	$(COMPILER) $(LINKER_OPTIONS) $(FUSSY_OPTIONS) $(SANITIZER_OPTIONS) -o $(EXECUTABLE) $(OBJECT_FILE) $(SDW_LINKER_FLAGS) $(SDL_LINKER_FLAGS)
	./$(EXECUTABLE)

# 高性能模式（用于渲染动画）
speedy: $(SDW_OBJECT_FILES)
	@mkdir -p $(BUILD_DIR)
	$(COMPILER) $(COMPILER_OPTIONS) $(SPEEDY_OPTIONS) -o $(OBJECT_FILE) $(SOURCE_FILE) $(SDL_COMPILER_FLAGS) $(SDW_COMPILER_FLAGS) $(GLM_COMPILER_FLAGS)
	$(COMPILER) $(LINKER_OPTIONS) $(SPEEDY_OPTIONS) -o $(EXECUTABLE) $(OBJECT_FILE) $(SDW_LINKER_FLAGS) $(SDL_LINKER_FLAGS)
	./$(EXECUTABLE)

# 生产模式
production: $(SDW_OBJECT_FILES)
	@mkdir -p $(BUILD_DIR)
	$(COMPILER) $(COMPILER_OPTIONS) -o $(OBJECT_FILE) $(SOURCE_FILE) $(SDL_COMPILER_FLAGS) $(SDW_COMPILER_FLAGS) $(GLM_COMPILER_FLAGS)
	$(COMPILER) $(LINKER_OPTIONS) -o $(EXECUTABLE) $(OBJECT_FILE) $(SDW_LINKER_FLAGS) $(SDL_LINKER_FLAGS)
	./$(EXECUTABLE)

# 编译 sdw 库
$(BUILD_DIR)/%.o: $(SDW_DIR)%.cpp
	@mkdir -p $(BUILD_DIR)
	$(COMPILER) $(COMPILER_OPTIONS) -c -o $@ $^ $(SDL_COMPILER_FLAGS) $(SDW_COMPILER_FLAGS) $(GLM_COMPILER_FLAGS)

# 清理
clean:
	rm -rf $(BUILD_DIR)/*
video:
	ffmpeg -framerate 60 -i frames/frame-%04d.ppm \
		-c:v libx264 -pix_fmt yuv420p output.mp4
	@echo "Video saved to output.mp4"

clean_frames:
	rm -f frames/*.ppm
	@echo "Frames cleaned."

# ==================== Cornell Box 50分版本 ====================
CORNELL_SOURCE := renderer-50.cpp
CORNELL_OBJECT := $(BUILD_DIR)/renderer-50.o
CORNELL_EXECUTABLE := $(BUILD_DIR)/renderer-50

cornell-debug: $(SDW_OBJECT_FILES)
	@mkdir -p $(BUILD_DIR)
	$(COMPILER) $(COMPILER_OPTIONS) $(DEBUG_OPTIONS) -o $(CORNELL_OBJECT) $(CORNELL_SOURCE) $(SDL_COMPILER_FLAGS) $(SDW_COMPILER_FLAGS) $(GLM_COMPILER_FLAGS)
	$(COMPILER) $(LINKER_OPTIONS) $(DEBUG_OPTIONS) -o $(CORNELL_EXECUTABLE) $(CORNELL_OBJECT) $(SDW_LINKER_FLAGS) $(SDL_LINKER_FLAGS)
	@echo "编译完成！运行: ./$(CORNELL_EXECUTABLE)"
	./$(CORNELL_EXECUTABLE)

cornell-speedy: $(SDW_OBJECT_FILES)
	@mkdir -p $(BUILD_DIR)
	$(COMPILER) $(COMPILER_OPTIONS) $(SPEEDY_OPTIONS) -o $(CORNELL_OBJECT) $(CORNELL_SOURCE) $(SDL_COMPILER_FLAGS) $(SDW_COMPILER_FLAGS) $(GLM_COMPILER_FLAGS)
	$(COMPILER) $(LINKER_OPTIONS) $(SPEEDY_OPTIONS) -o $(CORNELL_EXECUTABLE) $(CORNELL_OBJECT) $(SDW_LINKER_FLAGS) $(SDL_LINKER_FLAGS)
	@echo "编译完成！运行: ./$(CORNELL_EXECUTABLE)"
	./$(CORNELL_EXECUTABLE)
cornell-movie: clean_frames
	$(MAKE) cornell-speedy
	$(MAKE) video

.PHONY: default debug diagnostic speedy production clean cornell-debug cornell-speedy
