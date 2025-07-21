#include "Models.h"

void Models::registerModel(std::shared_ptr<BaseModel> model) {
    modelList.push_back(model);
}

std::shared_ptr<BaseModel> Models::getModel(const std::string& name) const {
    for (const auto& model : modelList) {
        if (model->getName() == name)
            return model;
    }
    return nullptr;
}

std::vector<std::string> Models::listModelNames() const {
    std::vector<std::string> names;
    for (const auto& model : modelList)
        names.push_back(model->getName());
    return names;
}
