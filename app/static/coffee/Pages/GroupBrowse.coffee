class GroupBrowse extends Page_Template
    model: () =>
        @data ?= new Data()
        @table ?= new Table(@data)
        return
    controller: (options) =>
        @model()
        return
    view: () =>
        return [
            super()
            @table.view()
        ]