#include <SFML/Graphics.hpp>
#include "Paint.h"
#include <thread>
#include <mutex>
#include <sstream>

sf::Image image;
std::mutex mutex;

void renderThread(const Scene & scene) {
    Camera camera([](int x, int y, unsigned char r, unsigned char g, unsigned char b) {
        mutex.lock();
        image.setPixel(x, y, {r, g, b, 255});
        mutex.unlock();
    }, image.getSize().x, image.getSize().y);
    camera.Render(scene);
}

int main() {
    sf::Texture texture;
    sf::RenderWindow window(sf::VideoMode(2560, 1600), "ProcessorRendering");
    window.setKeyRepeatEnabled(false);
    window.setActive();

    image.create(window.getSize().x, window.getSize().y);

    Mesh shipMesh("ship"), cubeMesh("cube"), dragonMesh("dragon"), girlMesh("girl");
    Camera camera([](int x, int y, unsigned char r, unsigned char g, unsigned char b) {
        image.setPixel(x, y, {r, g, b, 255});
    }, window.getSize().x, window.getSize().y);
    camera.setLightDirection({0, 0.5f, 1.f});
    Object ship1(shipMesh), cube1(cubeMesh), cube2(cubeMesh), dragon(dragonMesh), girl(girlMesh);
    ship1.transform.position = {0, -40, 250};
    ship1.transform.eulerAngles(-90, 0, 0);
    cube1.transform.position = {0, -1, 4};
    cube1.transform.eulerAngles(0, 0, 0);
    cube1.setColor(100, 200, 100);
    cube2.transform.position = {-1, -1, 4};
    cube2.transform.eulerAngles(0, 45, 45);
    cube2.setColor(200, 200, 0);
    dragon.transform.position = {50, -120, 500};
    dragon.transform.eulerAngles(0, 240, 0);
    dragon.setColor(200, 100, 50);
    girl.transform.position = {0, -40, 70};
    girl.transform.eulerAngles(0, 200, 0);

    Scene shipScene, cubeScene, dragonScene, girlScene;
    shipScene.addToScene(ship1);

    cubeScene.addToScene(cube1);
    cubeScene.addToScene(cube2);

    dragonScene.addToScene(dragon);

    girlScene.addToScene(girl);

    std::thread *t1;
    //std::thread t1(&renderThread, dragonScene);
    //t1.join();

    sf::Clock clock;
    int frameCount = 0;
    bool isKeyPressed = false, isRenderStarted = false;
    while (window.isOpen()) {
        ++frameCount;
        if (frameCount == 100) {
            float fps = 100/clock.restart().asSeconds();
            std::ostringstream sout;
            sout << fps;
            //window.setTitle(sout.str());
        }
        sf::Event event{};
        while (window.pollEvent(event)) {
            switch (event.type) {
                case sf::Event::Closed:
                    window.close();
                    t1->join();
                    return 0;
                case sf::Event::KeyPressed:
                    isKeyPressed = true;
                    break;
                default:
                    break;
            }
        }
        if (isKeyPressed) {
            //update ship scene
            ship1.transform.rotate({0.f, -1.f, 0.f}, false);
            //update cube scene
            cube1.transform.rotate({0.0f, 0.2f, 0.0f}, false);
            cube2.transform.rotate({0.2f, 0.2f, 0.0f}, false);
            cube2.transform.position += {0.005f, 0, 0};

            window.clear();
//            image.create(window.getSize().x, window.getSize().y);
//            camera.Render(shipScene);
//            image.flipVertically();
            if (!isRenderStarted) {
                isRenderStarted = true;
                t1 = new std::thread(&renderThread, dragonScene);
            }
            mutex.lock();
            sf::Image image1 = image;
            mutex.unlock();
            image1.flipVertically();
            texture.loadFromImage(image1);
            window.draw(sf::Sprite(texture));
            window.display();
        }
    }

    return 0;
}

