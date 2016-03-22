class UploadGenome extends PageTemplate
    model: () =>
        @form ?= new GenomeForm()
        return
    controller: (options) =>
        @model()
        return
    view: () =>
        return [
            super()
            @form.view()
        ]