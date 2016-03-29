class UploadGenome extends PageTemplate
    model: () =>
        @form = new UploadForm()
        return
    controller: (options) =>
        @model()
        return
    view: () =>
        return [
            super()
            @form.view()
        ]