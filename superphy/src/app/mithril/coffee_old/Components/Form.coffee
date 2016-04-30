class UploadForm
    model: () =>
        @genome_name = {'value': m.prop(''), 'label':'Genome Name'}
        @serotype = {'value': m.prop(''), 'label':'Serotype'}
        @isolation_host = {'value': m.prop(''), 'label':'Isolation Host'}
        @isolation_source = {'value': m.prop(''), 'label':'Isolation Source'}

    constructor:()->
        @model()

    send: () =>
        m.request(
            method: "POST",
            url: "http://#{location.hostname}:5000/example/echo"
            data: {
                genomes: [
                    'genome_name': @genome_name.value()
                    'serotype': @serotype.value()
                    'isolation_host': @isolation_host.value()
                    'isolation_source': @isolation_source.value()
                ]
            }
            datatype: 'json',
            type: (response)=>
                console.log(response)
        )
        console.log([
            @genome_name.value()
            @serotype.value()
            @isolation_host.value()
            @isolation_source.value()
        ])
    view: () =>
        return [


            m('div', @genome_name.label
                m('input[type=text]', {
                    oninput: m.withAttr('value', @genome_name.value)
                    value: @genome_name.value()
                }))
            m('div', @serotype.label
                m('input[type=text]', {
                    oninput: m.withAttr('value', @serotype.value)
                    value: @serotype.value()
                }))
            m('div', @isolation_host.label
                m('input[type=text]', {
                    oninput: m.withAttr('value', @isolation_host.value)
                    value: @isolation_host.value()
                }))
            m('div', @isolation_source.label
                m('input[type=text]', {
                    oninput: m.withAttr('value', @isolation_source.value)
                    value: @isolation_source.value()
                }))

            m("button[type=button]", {onclick: @send}, "Done")
        ]




class Uploader
    constructor: ()->
        @controller()
    upload: (options) ->
        formData = new formData
        for key in options.data
            for i in [0 .. options.data[key].length]
                formData.append(key, options.data[key][i])
        options.serialize (value) ->
            return value

        options.data = formData

        return m.request(options)
    serialize: (files) ->
        promises = files.map((file) ->
            deferred = m.deferred()
            reader = new FileReader
            reader.readAsDataURL()
            reader.onloadend = (e) ->
                deffered.resolve(e.result)
            reader.onerror = deffered.reject
            return deferered.promises
        )
        return m.sync(promises)
    controller: (args) ->
        @noop = (e) ->
            e.preventDefault()
        @update = (e) ->
            e.preventDefault()
            if (typeof args.onchange == "function")
                args.onchange([].slice.call((e.dataTransfer || e.target).files))
    view: ()=>
        return m(".uploader", {
            ondragover: @controller.noop
            ondrop: @controller.update
            })

"""
class Foo
    @controller: (args) ->
        @bang = (e) ->
            console.log("Bang")
        return
    @view: (ctrl) ->
        return m("button", {onclick: ctrl.bang}, "Bang")

Test = {
    controller: (args) ->
        @bing = (e) ->
            console.log("bing")
        return
    view: (ctrl) =>
        return m("button", {onclick: ctrl.bing}, "Bing")
}
"""