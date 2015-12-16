class GroupBrowse extends Page_Template
    model: () =>
        @search = m.prop('')
        @value = m.prop('')
        @data ?= new Data()
        @table ?= new Table(@data)
        return
    controller: (options) =>
        @model()
        return
    view: () =>
        return [
            super()
            m('input', {oninput: m.withAttr('value', @value), value : @value() })
            m('button', { disabled : !@value(), onclick: => if @value() then @search(@value())},["Search"])
            @table.view(@search())
        ]