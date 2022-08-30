
#include "Engine.hpp"
#include "MatrixSizeExecutor.hpp"
#include "UpdateCountExecutor.hpp"

namespace {
  std::unique_ptr<EngineModeExecutor> makeExecutor(EngineMode mode, Engine *engine) {
    // Create the correct executor for the given mode
    // Note that the executor can use the engine to get access to the command line arguments
    switch (mode) {
      case EngineMode::kMatrixSize:
        return std::make_unique<MatrixSizeExecutor>(engine);
      case EngineMode::kUpdateCount:
        return std::make_unique<UpdateCountExecutor>(engine);
      default:
        throw std::runtime_error("Unknown engine mode");
    }
  }

  void checkOutputPath(std::filesystem::path output_path) {

    if (not std::filesystem::exists(output_path)) {
      output_path = std::filesystem::absolute(output_path);
      // ensure that the parent path exists
      std::filesystem::create_directories(output_path.parent_path());
      return;
    }

    char answer = '\0';
    while (answer != 'y' and answer != 'n') {
      std::cout << "Warning ! A file named " << output_path
                << " Already exists !\n\tOverwrite ? (y/n) ";
      std::string buffer;
      std::cin >> buffer;
      answer = buffer[0];
    }

    if (answer == 'n') { throw std::runtime_error("Aborted"); }

    std::filesystem::remove(output_path);
  }
}// namespace

Engine::Engine(int argc, char **argv) : argv(argv), argc(argc) {
  if (argc < 3) {
    std::cout << "Missing arguments" << std::endl;
    std::cout << "Usage: " << argv[0] << " <dataset_output_path> "
              << "<update|matrix_size|...> (<args>)" << std::endl;
    std::abort();
  }

  if (std::string(argv[2]) == "update") {
    mode = EngineMode::kUpdateCount;
  } else if (std::string(argv[2]) == "matrix_size") {
    mode = EngineMode::kMatrixSize;
  } else {
    std::cout << "Error: Unknown mode: " << argv[2] << std::endl;
    std::abort();
  }

  dataset_output_path = argv[1];
  checkOutputPath(dataset_output_path);
}

int Engine::exec() {
  auto executor = makeExecutor(mode, this);
  auto return_code = executor->exec();

  return return_code;
}
