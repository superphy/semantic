/**
 * Created by Stephen Kan on 20/10/15.
 */

describe("Test hello.js", function() {
    it("return 'hello'", function () {
        var result = Hello.world();
        expect(result).toBe("Hello");
    });
});